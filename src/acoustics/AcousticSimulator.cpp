#include "acoustics/AcousticSimulator.hpp"

#include <algorithm>
#include <numbers>

namespace mspatial::acoustics {
namespace {

constexpr float kValidationTolerance = 5e-3f;
constexpr float kPi = std::numbers::pi_v<float>;
constexpr float kGolden = 0.6180339887f;
constexpr std::uint32_t kMaxPortalChain = 8;

AcousticRayRecord segment(Vec3 from, Vec3 direction, float tMin, float tMax) {
    return {from.x, from.y, from.z, tMin, direction.x, direction.y, direction.z, tMax};
}

Vec3 fibonacciSphere(std::uint32_t index, std::uint32_t count) {
    const float z = 1.f - (2.f * float(index) + 1.f) / float(count);
    const float r = std::sqrt(std::max(0.f, 1.f - z * z));
    const float phi = 2.f * kPi * std::fmod(float(index) * kGolden, 1.f);
    return vec3(r * std::cos(phi), r * std::sin(phi), z);
}

float energy(const Bands& amplitude) {
    float sum = 0.f;
    for (float a : amplitude) sum += a * a;
    return sum;
}

float peak(const Bands& amplitude) { return *std::max_element(amplitude.begin(), amplitude.end()); }

// R2 low-discrepancy point in the unit square.
void r2(std::uint32_t index, float& u, float& v) {
    // Fixed point, so the sequence keeps its quality for any index.
    u = float(index * 3242174889u) / 4294967296.f;
    v = float(index * 2447445414u) / 4294967296.f;
}

bool samePlane(const WorldPlane& a, const WorldPlane& b) {
    const float d = dot(a.normal, b.normal);
    return (d > 0.9995f && std::fabs(a.offset - b.offset) < 2e-3f) || (d < -0.9995f && std::fabs(a.offset + b.offset) < 2e-3f);
}

}  // namespace

std::uint64_t PathKey::hash() const {
    std::uint64_t h = 0xcbf29ce484222325ull ^ order;
    for (std::uint32_t i = 0; i < order; ++i) {
        for (std::uint32_t value : {facets[i].instance, facets[i].facet}) {
            h = (h ^ value) * 0x100000001b3ull;
            h ^= h >> 29;
        }
    }
    return h | 1ull;  // Zero marks an empty slot.
}

bool solveImagePath(const Scene& scene, const PathKey& key, Vec3 source, Vec3 listener, const Bands& air, ImagePath& out) {
    if (key.order == 0 || key.order > kMaxReflectionOrder) return false;
    std::array<WorldPlane, kMaxReflectionOrder> planes;
    std::array<Vec3, kMaxReflectionOrder> images;
    Vec3 image = source;
    for (std::uint32_t k = 0; k < key.order; ++k) {
        if (!scene.facetValid(key.facets[k].instance, key.facets[k].facet)) return false;
        planes[k] = scene.worldPlane(key.facets[k].instance, key.facets[k].facet);
        image = image - planes[k].normal * (2.f * (dot(planes[k].normal, image) - planes[k].offset));
        images[k] = image;
    }
    out.length = length(listener - image);
    if (out.length < 1e-4f) return false;
    Vec3 target = listener;
    for (std::uint32_t k = key.order; k-- > 0;) {
        const float targetSide = dot(planes[k].normal, target) - planes[k].offset;
        const float imageSide = dot(planes[k].normal, images[k]) - planes[k].offset;
        if (targetSide * imageSide >= 0.f) return false;  // A real reflection crosses the plane on its way to the image.
        const float s = targetSide / (targetSide - imageSide);
        out.points[k] = target + (images[k] - target) * s;
        target = out.points[k];
    }
    // Every side test passed, so each image lay across its plane from the next point: that also proves the previous image, and at the end
    // the source, is on the listener's side of each reflection.
    out.emission = normalize(out.points[0] - source);
    out.arrival = normalize(out.points[key.order - 1] - listener);
    for (std::size_t band = 0; band < kBandCount; ++band) {
        float reflected = 1.f;
        for (std::uint32_t k = 0; k < key.order; ++k) {
            const Material& material = scene.facetMaterial(key.facets[k].instance, key.facets[k].facet);
            reflected *= (1.f - material.absorption[band]) * (1.f - material.scattering[band]);
        }
        out.amplitude[band] = std::sqrt(reflected * std::exp(-air[band] * out.length)) / std::max(out.length, 0.1f);
    }
    return true;
}

Simulator::Simulator(Scene& scene, TraceBackend& backend, const SimulatorConfig& config) : scene_(scene), backend_(backend), config_(config) {
    config_.directRays = std::max<std::uint32_t>(config_.directRays, 1);
    config_.sequenceCacheSize = std::max<std::uint32_t>(config_.sequenceCacheSize, 64);
    sequences_.resize(config_.sequenceCacheSize);
    sequenceIndex_.reserve(config_.sequenceCacheSize * 2);
}

std::uint32_t Simulator::addSource(const SourceDesc& desc) {
    for (std::uint32_t i = 0; i < kMaxSources; ++i) {
        if (sources_[i].active) continue;
        SourceState& s = sources_[i];
        const std::uint32_t generation = s.generation + 1;
        s = SourceState{};
        s.generation = generation;
        s.desc = desc;
        s.active = true;
        return i;
    }
    return kNone;
}

void Simulator::setSource(std::uint32_t id, const SourceDesc& desc) { sources_.at(id).desc = desc; }

void Simulator::removeSource(std::uint32_t id) {
    sources_.at(id).active = false;
    ++sources_[id].generation;
}

void Simulator::tick(float dt) {
    ++tick_;
    backend_.syncScene(scene_);
    roomStates_.resize(scene_.rooms().size());
    portalStates_.resize(scene_.portals().size());
    while (TraceBatch* batch = backend_.pollCompleted()) {
        consume(*batch);
        backend_.release(batch);
    }
    if (TraceBatch* batch = backend_.acquireBatch()) {
        prepare(*batch);
        backend_.submit(batch);
    }
    publish(dt);
}

bool Simulator::pushPathJob(TraceBatch& batch, Jobs& jobs, std::uint32_t source, const PathKey& key, const ImagePath& path, bool candidate) {
    const std::uint32_t rays = key.order + 1;
    if (batch.closestCount + rays > backend_.capacity().closestRays) return false;
    const SourceState& state = sources_[source];
    PathJob job{source, state.generation, batch.closestCount, key, {}, candidate};
    const float eps = backend_.selfHitOffset();
    Vec3 from = state.desc.position;
    for (std::uint32_t i = 0; i < rays; ++i) {
        const Vec3 to = i < key.order ? path.points[i] : listener_;
        const float distance = length(to - from);
        job.lengths[i] = distance;
        const float tMax = i < key.order ? distance + kValidationTolerance : std::max(distance - kValidationTolerance, eps);
        batch.closestRays[batch.closestCount++] = segment(from, (to - from) / std::max(distance, 1e-6f), eps, tMax);
        from = to;
    }
    jobs.paths.push_back(job);
    return true;
}

void Simulator::prepare(TraceBatch& batch) {
    Jobs& jobs = jobs_[batch.slot];
    jobs.paths.clear();
    jobs.direct.clear();
    jobs.portals.clear();
    jobs.probes.clear();
    jobs.probeDirections.clear();
    const TraceCapacity& capacity = backend_.capacity();
    const float eps = backend_.selfHitOffset();

    // Portals and rooms come first: they are few, and they steer the whole late field.
    const auto& portals = scene_.portals();
    for (std::uint32_t p = 0; p < portals.size(); ++p) {
        PortalState& state = portalStates_[p];
        const PortalDesc& portal = portals[p];
        for (std::uint32_t k = 0; k < config_.portalRaysPerTick && batch.closestCount < capacity.closestRays; ++k) {
            const std::uint32_t sample = state.cursor++ % state.samples.size();
            float u, v;
            r2(sample, u, v);
            const Vec3 point = portal.center + portal.halfU * (2.f * u - 1.f) + portal.halfV * (2.f * v - 1.f);
            jobs.portals.push_back({p, sample, batch.closestCount});
            batch.closestRays[batch.closestCount++] = segment(point - portal.normal * portal.halfThickness, portal.normal, 0.f, 2.f * portal.halfThickness);
        }
    }
    const auto& rooms = scene_.rooms();
    for (std::uint32_t r = 0; r < rooms.size(); ++r) {
        if (rooms[r].volume > 0.f && rooms[r].surfaceArea > 0.f) continue;
        RoomState& state = roomStates_[r];
        state.samples.resize(config_.probeWindow);
        for (std::uint32_t k = 0; k < config_.probeRaysPerTick && batch.closestCount < capacity.closestRays; ++k) {
            const std::uint32_t sample = state.cursor++ % config_.probeWindow;
            const Vec3 direction = fibonacciSphere(sample, config_.probeWindow);
            jobs.probes.push_back({r, sample, batch.closestCount});
            jobs.probeDirections.push_back(direction);
            batch.closestRays[batch.closestCount++] = segment(rooms[r].probe, direction, 0.f, config_.probeMaxDistance);
        }
    }

    // Direct line: a few layered rays spread over the source disc, the pattern shifting each tick so smoothing averages more of the disc.
    for (std::uint32_t s = 0; s < kMaxSources; ++s) {
        SourceState& source = sources_[s];
        if (!source.active) continue;
        const Vec3 toSource = source.desc.position - listener_;
        const float distance = length(toSource);
        if (distance < 1e-4f || batch.layeredCount + config_.directRays > capacity.layeredRays) continue;
        const Vec3 axis = toSource / distance;
        const Vec3 helper = std::fabs(axis.z) < 0.9f ? vec3(0, 0, 1) : vec3(1, 0, 0);
        const Vec3 tangent = normalize(cross(axis, helper));
        const Vec3 bitangent = cross(axis, tangent);
        jobs.direct.push_back({s, source.generation, batch.layeredCount, config_.directRays});
        for (std::uint32_t k = 0; k < config_.directRays; ++k) {
            Vec3 target = source.desc.position;
            if (k > 0) {
                float u, v;
                r2(std::uint32_t(tick_) * config_.directRays + k, u, v);
                const float radius = source.desc.radius * std::sqrt(u), angle = 2.f * kPi * v;
                target = target + tangent * (radius * std::cos(angle)) + bitangent * (radius * std::sin(angle));
            }
            const Vec3 line = target - listener_;
            const float span = length(line);
            batch.layeredRays[batch.layeredCount++] = segment(listener_, line / span, eps, span);
        }
    }

    // Early paths share what is left. Each source gets an equal quota, and any it leaves unused passes to the next: the start rotates, so a tight
    // budget stays fair. Inside a quota, confirming kept paths comes first but may not use all of it, or a full reservoir could never improve.
    std::uint32_t sourcesLeft = 0;
    for (const SourceState& source : sources_) sourcesLeft += source.active;
    fairSource_ = (fairSource_ + 1) % kMaxSources;
    for (std::uint32_t i = 0; i < kMaxSources; ++i) {
        const std::uint32_t s = (fairSource_ + i) % kMaxSources;
        SourceState& source = sources_[s];
        if (!source.active) continue;
        const std::uint32_t quota = (capacity.closestRays - batch.closestCount) / sourcesLeft--;
        const std::uint32_t quotaEnd = batch.closestCount + quota;

        // Geometry is exact for every kept path every tick. Rays only confirm visibility, so a path skipped this tick still tracks the live scene.
        std::uint32_t freeSlots = 0;
        float weakest = INFINITY;
        std::array<ImagePath, kMaxPathsPerSource> solved;
        std::array<bool, kMaxPathsPerSource> solvable{};
        for (std::uint32_t k = 0; k < kMaxPathsPerSource; ++k) {
            PathSlot& slot = source.slots[k];
            if (slot.id == 0) {
                ++freeSlots;
                continue;
            }
            solvable[k] = solveImagePath(scene_, slot.key, source.desc.position, listener_, config_.air, solved[k]);
            if (!solvable[k]) {
                slot.visible = false;
                continue;
            }
            slot.path = solved[k];
            if (slot.visible) weakest = std::min(weakest, energy(slot.path.amplitude));
        }
        const std::uint32_t confirmEnd = sequenceCount_ > 0 ? batch.closestCount + quota - quota / 4 : quotaEnd;
        for (std::uint32_t k = 0; k < kMaxPathsPerSource; ++k) {
            const std::uint32_t index = (source.confirmCursor + k) % kMaxPathsPerSource;
            if (!solvable[index] || source.slots[index].retiring) continue;
            if (batch.closestCount + source.slots[index].key.order + 1 > confirmEnd) {
                source.confirmCursor = index;
                break;
            }
            pushPathJob(batch, jobs, s, source.slots[index].key, solved[index], false);
        }

        // New candidates: the analytic filter is cheap, so each source reads a slice of the shared sequence cache and validates the strongest it saw.
        if (sequenceCount_ == 0) continue;
        struct Candidate {
            const PathKey* key;
            ImagePath path;
            float strength;
        };
        std::array<Candidate, 16> best;
        const std::uint32_t keep = std::min<std::uint32_t>(config_.newValidationsPerSource, std::uint32_t(best.size()));
        std::uint32_t held = 0;
        const std::uint32_t checks = std::min<std::uint32_t>(config_.sequenceChecksPerSource, std::uint32_t(sequenceCount_));
        for (std::uint32_t c = 0; c < checks; ++c) {
            const PathKey& key = sequences_[source.cursor++ % sequenceCount_];
            const std::uint64_t id = key.hash();
            bool known = false;
            for (const PathSlot& slot : source.slots) known = known || slot.id == id;
            if (known) continue;
            ImagePath path;
            if (!solveImagePath(scene_, key, source.desc.position, listener_, config_.air, path)) continue;
            if (peak(path.amplitude) < config_.minimumPathGain) continue;
            const float strength = energy(path.amplitude);
            if (freeSlots == 0 && strength < weakest * config_.replacementMargin) continue;
            if (held == keep && strength <= best[held - 1].strength) continue;
            bool onFacets = true;
            for (std::uint32_t k = 0; k < key.order && onFacets; ++k)
                onFacets = scene_.nearFacet(key.facets[k].instance, key.facets[k].facet, path.points[k], 0.02f);
            if (!onFacets) continue;
            std::uint32_t at = std::min(held, keep - 1);
            while (at > 0 && best[at - 1].strength < strength) {
                best[at] = best[at - 1];
                --at;
            }
            best[at] = {&key, path, strength};
            held = std::min(held + 1, keep);
        }
        for (std::uint32_t c = 0; c < held; ++c) {
            if (batch.closestCount + best[c].key->order + 1 > quotaEnd) break;
            pushPathJob(batch, jobs, s, *best[c].key, best[c].path, true);
        }
    }
    issuedClosest_ = batch.closestCount;

    jobs.discovery = std::min(config_.discoveryRays, capacity.discoveryRays);
    for (std::uint32_t i = 0; i < jobs.discovery; ++i) {
        // One R2 sequence runs across ticks, so the whole sphere is covered evenly however few rays each tick has.
        float u, v;
        r2(std::uint32_t(tick_ * jobs.discovery) + i, u, v);
        const float z = 1.f - 2.f * u, ring = std::sqrt(std::max(0.f, 1.f - z * z)), phi = 2.f * kPi * v;
        const Vec3 d = vec3(ring * std::cos(phi), ring * std::sin(phi), z);
        batch.discoverySeeds[i] = {listener_.x, listener_.y, listener_.z, std::uint32_t(tick_) * 2654435761u + i, d.x, d.y, d.z, 0};
    }
    batch.discoveryCount = jobs.discovery;
}

void Simulator::consume(const TraceBatch& batch) {
    if (batch.failed) return;
    const Jobs& jobs = jobs_[batch.slot];
    consumeDiscovery(batch, jobs);
    consumePaths(batch, jobs);
    consumeDirect(batch, jobs);
    consumeSamples(batch, jobs);
}

void Simulator::rememberSequence(const PathKey& key) {
    const std::uint64_t id = key.hash();
    if (sequenceIndex_.contains(id)) return;
    if (sequenceCount_ == sequences_.size()) sequenceIndex_.erase(sequences_[sequenceHead_].hash());
    sequences_[sequenceHead_] = key;
    sequenceIndex_[id] = sequenceHead_;
    sequenceHead_ = (sequenceHead_ + 1) % std::uint32_t(sequences_.size());
    sequenceCount_ = std::min(sequenceCount_ + 1, sequences_.size());
}

void Simulator::consumeDiscovery(const TraceBatch& batch, const Jobs& jobs) {
    const std::uint32_t maxOrder = std::min(backend_.capacity().maxOrder, kMaxReflectionOrder);
    for (std::uint32_t i = 0; i < batch.discoveryCount && i < jobs.discovery; ++i) {
        const AcousticDiscoveryVertex* chain = batch.discoveryVertices + i * backend_.capacity().maxOrder;
        std::array<FacetRef, kMaxReflectionOrder> walked;
        std::uint32_t count = 0;
        for (; count < maxOrder && chain[count].instance != ACOUSTIC_MISS; ++count) {
            walked[count] = {chain[count].instance, chain[count].facet};
            if (count > 0 && walked[count] == walked[count - 1]) break;  // Two reflections off one plane cannot happen.
            // Every prefix is a sequence, ending at the listener's first surface. The source side comes first in the key.
            PathKey key;
            key.order = count + 1;
            for (std::uint32_t k = 0; k <= count; ++k) key.facets[k] = walked[count - k];
            rememberSequence(key);
        }
    }
}

void Simulator::admit(SourceState& source, const PathKey& key, const ImagePath& path) {
    const std::uint64_t id = key.hash();
    PathSlot* free = nullptr;
    PathSlot* weakest = nullptr;
    float weakestEnergy = INFINITY;
    for (PathSlot& slot : source.slots) {
        if (slot.id == id) return;
        if (slot.id == 0) {
            if (!free) free = &slot;
        } else if (slot.visible && energy(slot.path.amplitude) < weakestEnergy) {
            weakestEnergy = energy(slot.path.amplitude);
            weakest = &slot;
        }
    }
    if (!free) {
        // The weakest path fades out first; the newcomer takes its slot on a later pass, so no slot ever changes path at full level.
        if (weakest && energy(path.amplitude) > weakestEnergy * config_.replacementMargin) {
            weakest->visible = false;
            weakest->retiring = true;
        }
        return;
    }
    *free = PathSlot{};
    free->key = key;
    free->id = id;
    free->path = path;
    free->visible = true;
}

void Simulator::consumePaths(const TraceBatch& batch, const Jobs& jobs) {
    for (const PathJob& job : jobs.paths) {
        SourceState& source = sources_[job.source];
        if (!source.active || source.generation != job.generation) continue;
        bool clear = true;
        for (std::uint32_t i = 0; i <= job.key.order && clear; ++i) {
            const AcousticHitRecord& hit = batch.closestHits[job.firstRay + i];
            if (i == job.key.order) {
                clear = hit.instance == ACOUSTIC_MISS;
                break;
            }
            const FacetRef expected = job.key.facets[i];
            if (hit.instance == ACOUSTIC_MISS || hit.t < job.lengths[i] - kValidationTolerance || !scene_.facetValid(expected.instance, expected.facet)) {
                clear = false;
            } else if (hit.instance != expected.instance || hit.facet != expected.facet) {
                // Modular walls: a neighbouring instance in the same plane is the same reflector.
                clear = scene_.facetValid(hit.instance, hit.facet) &&
                        samePlane(scene_.worldPlane(hit.instance, hit.facet), scene_.worldPlane(expected.instance, expected.facet));
            }
        }
        const std::uint64_t id = job.key.hash();
        if (job.candidate) {
            if (clear) {
                ImagePath path;
                if (solveImagePath(scene_, job.key, source.desc.position, listener_, config_.air, path)) admit(source, job.key, path);
            }
            continue;
        }
        for (PathSlot& slot : source.slots) {
            if (slot.id != id) continue;
            if (clear) {
                slot.misses = 0;
                slot.visible = !slot.retiring;
            } else if (++slot.misses >= config_.missTolerance) {
                slot.visible = false;
            }
            break;
        }
    }
}

void Simulator::consumeDirect(const TraceBatch& batch, const Jobs& jobs) {
    for (const DirectJob& job : jobs.direct) {
        SourceState& source = sources_[job.source];
        if (!source.active || source.generation != job.generation) continue;
        Bands sum = uniformBands(0.f);
        float clear = 0.f;
        for (std::uint32_t r = 0; r < job.rays; ++r) {
            const AcousticLayeredHit& layers = batch.layeredHits[job.firstRay + r];
            if (layers.count == 0) clear += 1.f;
            Bands through = uniformBands(1.f);
            for (std::uint32_t k = 0; k < layers.count; ++k) {
                if (!scene_.facetValid(layers.hits[k].instance, layers.hits[k].facet)) continue;
                const Material& material = scene_.facetMaterial(layers.hits[k].instance, layers.hits[k].facet);
                for (std::size_t b = 0; b < kBandCount; ++b) through[b] *= material.transmission[b];
            }
            for (std::size_t b = 0; b < kBandCount; ++b) sum[b] += through[b];
        }
        for (std::size_t b = 0; b < kBandCount; ++b) source.directTarget[b] = sum[b] / float(job.rays);
        source.visibilityTarget = clear / float(job.rays);
    }
}

void Simulator::consumeSamples(const TraceBatch& batch, const Jobs& jobs) {
    for (const SampleJob& job : jobs.portals) {
        if (job.owner >= portalStates_.size()) continue;
        PortalState& state = portalStates_[job.owner];
        const AcousticHitRecord& hit = batch.closestHits[job.ray];
        Bands open = uniformBands(1.f);
        if (hit.instance != ACOUSTIC_MISS && scene_.facetValid(hit.instance, hit.facet)) open = scene_.facetMaterial(hit.instance, hit.facet).transmission;
        state.samples[job.sample] = open;
        state.filled[job.sample] = true;
    }
    for (PortalState& state : portalStates_) {
        Bands sum = uniformBands(0.f);
        float count = 0.f;
        for (std::size_t i = 0; i < state.samples.size(); ++i) {
            if (!state.filled[i]) continue;
            count += 1.f;
            for (std::size_t b = 0; b < kBandCount; ++b) sum[b] += state.samples[i][b];
        }
        if (count > 0.f)
            for (std::size_t b = 0; b < kBandCount; ++b) state.open[b] = sum[b] / count;
    }

    const auto& portals = scene_.portals();
    for (std::size_t i = 0; i < jobs.probes.size(); ++i) {
        const SampleJob& job = jobs.probes[i];
        if (job.owner >= roomStates_.size() || job.sample >= roomStates_[job.owner].samples.size()) continue;
        const AcousticHitRecord& hit = batch.closestHits[job.ray];
        const Vec3 direction = jobs.probeDirections[i];
        const Vec3 origin = scene_.rooms()[job.owner].probe;
        float range = hit.instance == ACOUSTIC_MISS ? config_.probeMaxDistance : hit.t;
        // A ray leaving through a portal stops counting at the opening: that area belongs to the coupling, not to the walls.
        bool leftThroughPortal = false;
        for (const PortalDesc& portal : portals) {
            if (portal.roomA != job.owner && portal.roomB != job.owner) continue;
            const float along = dot(direction, portal.normal);
            if (std::fabs(along) < 1e-5f) continue;
            const float t = dot(portal.center - origin, portal.normal) / along;
            if (t <= 0.f || t >= range) continue;
            const Vec3 local = origin + direction * t - portal.center;
            const float uu = dot(portal.halfU, portal.halfU), vv = dot(portal.halfV, portal.halfV);
            if (std::fabs(dot(local, portal.halfU)) <= uu && std::fabs(dot(local, portal.halfV)) <= vv) {
                range = t;
                leftThroughPortal = true;
            }
        }
        ProbeSample sample;
        sample.filled = true;
        sample.cubedRange = range * range * range;
        const float solidAngle = 4.f * kPi / float(config_.probeWindow);
        if (leftThroughPortal) {
            // Nothing to add: open area is carried by the link.
        } else if (hit.instance == ACOUSTIC_MISS || !scene_.facetValid(hit.instance, hit.facet)) {
            sample.area = range * range * solidAngle;  // Open to the outside: energy leaves as through a window.
            sample.absorption = uniformBands(sample.area);
        } else {
            const WorldPlane plane = scene_.worldPlane(hit.instance, hit.facet);
            sample.area = range * range * solidAngle / std::max(std::fabs(dot(direction, plane.normal)), 0.1f);
            const Material& material = scene_.facetMaterial(hit.instance, hit.facet);
            for (std::size_t b = 0; b < kBandCount; ++b) sample.absorption[b] = sample.area * material.absorption[b];
        }
        roomStates_[job.owner].samples[job.sample] = sample;
    }
}

void Simulator::updateRooms(AcousticSnapshot& out) {
    const auto& rooms = scene_.rooms();
    const auto& portals = scene_.portals();
    out.roomCount = std::min<std::uint32_t>(std::uint32_t(rooms.size()), kMaxRooms);
    out.linkCount = 0;
    for (std::uint32_t r = 0; r < out.roomCount; ++r) {
        RoomAcoustics& room = out.rooms[r];
        room = RoomAcoustics{};
        float volume = rooms[r].volume, surface = rooms[r].surfaceArea;
        Bands absorption = rooms[r].absorptionArea;
        if (volume <= 0.f || surface <= 0.f) {
            double cubed = 0.0, area = 0.0;
            std::array<double, kBandCount> measured{};
            std::uint32_t filled = 0;
            for (const ProbeSample& sample : roomStates_[r].samples) {
                if (!sample.filled) continue;
                ++filled;
                cubed += sample.cubedRange;
                area += sample.area;
                for (std::size_t b = 0; b < kBandCount; ++b) measured[b] += sample.absorption[b];
            }
            if (filled == 0) continue;  // Nothing known yet; the room stays inactive one more tick.
            const double scale = double(config_.probeWindow) / filled;
            if (volume <= 0.f) volume = float(cubed / filled * 4.0 * double(kPi) / 3.0);
            if (surface <= 0.f) {
                surface = float(area * scale);
                for (std::size_t b = 0; b < kBandCount; ++b) absorption[b] = float(measured[b] * scale);
            }
        }
        room.active = volume > 0.f;
        room.volume = volume;
        room.surfaceArea = surface;
        room.absorptionArea = absorption;
        room.meanFreePath = surface > 0.f ? 4.f * volume / surface : 0.f;
    }
    Bands opening[kMaxRooms];
    for (std::uint32_t r = 0; r < out.roomCount; ++r) opening[r] = uniformBands(0.f);
    for (std::uint32_t p = 0; p < portals.size(); ++p) {
        const PortalDesc& portal = portals[p];
        if (portal.roomA >= out.roomCount || portal.roomB >= out.roomCount) continue;
        if (!out.rooms[portal.roomA].active || !out.rooms[portal.roomB].active || out.linkCount + 2 > kMaxRoomLinks) continue;
        const float area = portal.area();
        const float forward = config_.speedOfSound * area / (4.f * out.rooms[portal.roomA].volume);
        const float backward = config_.speedOfSound * area / (4.f * out.rooms[portal.roomB].volume);
        RoomLink& ab = out.links[out.linkCount++];
        RoomLink& ba = out.links[out.linkCount++];
        ab.from = portal.roomA;
        ab.to = portal.roomB;
        ba.from = portal.roomB;
        ba.to = portal.roomA;
        for (std::size_t b = 0; b < kBandCount; ++b) {
            const float open = portalStates_[p].open[b];
            ab.rate[b] = forward * open;
            ba.rate[b] = backward * open;
            opening[portal.roomA][b] += area * open;
            opening[portal.roomB][b] += area * open;
        }
    }
    for (std::uint32_t r = 0; r < out.roomCount; ++r) {
        RoomAcoustics& room = out.rooms[r];
        if (!room.active) continue;
        for (std::size_t b = 0; b < kBandCount; ++b) {
            const float equivalent = equivalentAbsorption(config_.formula, room.surfaceArea, room.absorptionArea[b]);
            room.decayRate[b] = decayRate(room.volume, equivalent + opening[r][b], config_.air[b], config_.speedOfSound);
            room.t60[b] = t60FromRate(room.decayRate[b]);
        }
    }
}

void Simulator::updatePortalPaths(std::uint32_t sourceIndex, SourceAcoustics& out) {
    for (PortalPath& path : out.portalPaths) path = PortalPath{};
    const auto& portals = scene_.portals();
    const std::uint32_t room = out.room, listenerRoom = scene_.roomAt(listener_);
    if (room == kNone || listenerRoom == kNone || room == listenerRoom) return;
    const Vec3 sourcePosition = sources_[sourceIndex].desc.position;

    // Directed portals d: 2p runs a to b, 2p+1 runs b to a. portalDistance_ was filled from the listener in publish().
    struct Candidate {
        float distance;
        std::uint32_t directed;
    };
    std::array<Candidate, kMaxPortalPaths> best;
    std::uint32_t found = 0;
    for (std::uint32_t d = 0; d < portalDistance_.size(); ++d) {
        const PortalDesc& portal = portals[d / 2];
        if ((d & 1 ? portal.roomA : portal.roomB) != room || portalDistance_[d] == INFINITY) continue;
        const Candidate c{portalDistance_[d] + length(portal.center - sourcePosition), d};
        std::uint32_t at = std::min(found, kMaxPortalPaths);
        while (at > 0 && best[at - 1].distance > c.distance) --at;
        if (at >= kMaxPortalPaths) continue;
        for (std::uint32_t i = std::min(found, kMaxPortalPaths - 1); i > at; --i) best[i] = best[i - 1];
        best[at] = c;
        found = std::min(found + 1, kMaxPortalPaths);
    }
    for (std::uint32_t entry = 0; entry < found; ++entry) {
        // Chain in source-to-listener order. Following previous pointers from the source side walks toward the listener.
        std::array<std::uint32_t, kMaxPortalChain> chain;
        std::uint32_t count = 0;
        for (std::uint32_t d = best[entry].directed; d != kNone && count < kMaxPortalChain; d = portalPrevious_[d]) chain[count++] = d / 2;

        // String pulling: each crossing point moves toward the straight line between its neighbours, inside the opening.
        std::array<Vec3, kMaxPortalChain + 2> points;
        points[0] = sourcePosition;
        for (std::uint32_t i = 0; i < count; ++i) points[i + 1] = portals[chain[i]].center;
        points[count + 1] = listener_;
        for (int sweep = 0; sweep < 3; ++sweep) {
            for (std::uint32_t i = 0; i < count; ++i) {
                const PortalDesc& portal = portals[chain[i]];
                const Vec3 a = points[i], b = points[i + 2];
                const float sa = dot(a - portal.center, portal.normal), sb = dot(b - portal.center, portal.normal);
                Vec3 local = sa * sb < 0.f ? a + (b - a) * (sa / (sa - sb)) - portal.center : points[i + 1] - portal.center;
                const float uu = dot(portal.halfU, portal.halfU), vv = dot(portal.halfV, portal.halfV);
                const float u = std::clamp(dot(local, portal.halfU) / uu, -1.f, 1.f), v = std::clamp(dot(local, portal.halfV) / vv, -1.f, 1.f);
                points[i + 1] = portal.center + portal.halfU * u + portal.halfV * v;
            }
        }
        std::array<ApertureLeg, kMaxPortalChain> legs;
        float total = 0.f;
        for (std::uint32_t i = 0; i < count; ++i) {
            const PortalDesc& portal = portals[chain[i]];
            const Vec3 in = points[i + 1] - points[i], outward = points[i + 2] - points[i + 1];
            const float d = length(in), e = length(outward);
            legs[i].distance = d;
            legs[i].cosIn = std::clamp(std::fabs(dot(in, portal.normal)) / std::max(d, 1e-4f), 0.02f, 1.f);
            legs[i].cosOut = std::clamp(std::fabs(dot(outward, portal.normal)) / std::max(e, 1e-4f), 0.02f, 1.f);
            legs[i].area = portal.area();
            legs[i].open = portalStates_[chain[i]].open;
            total += d;
        }
        const float lastDistance = length(points[count + 1] - points[count]);
        total += lastDistance;
        const Bands energyGain = apertureChainEnergy(std::span(legs.data(), count), lastDistance, config_.speedOfSound);
        PortalPath& path = out.portalPaths[entry];
        path.portals = count;
        path.delaySeconds = total / config_.speedOfSound;
        path.excessSeconds = std::max(total - out.distance, 0.f) / config_.speedOfSound;
        path.arrival = normalize(points[count] - listener_);
        path.emission = normalize(points[1] - sourcePosition);
        for (std::size_t b = 0; b < kBandCount; ++b) path.gain[b] = std::sqrt(energyGain[b] * std::exp(-config_.air[b] * total));
    }
}

void Simulator::publish(float dt) {
    AcousticSnapshot& out = exchange_.writeSlot();
    out.tick = tick_;
    out.listener = listener_;
    out.listenerRoom = scene_.roomAt(listener_);
    updateRooms(out);

    // Shortest open route from the listener to every directed portal. Few portals, so the plain quadratic form is the right one.
    const auto& portals = scene_.portals();
    const std::size_t directed = portals.size() * 2;
    portalDistance_.assign(directed, INFINITY);
    portalPrevious_.assign(directed, kNone);
    portalDone_.assign(directed, 0);
    const auto entered = [&](std::uint32_t d) { return d & 1 ? portals[d / 2].roomA : portals[d / 2].roomB; };
    const auto left = [&](std::uint32_t d) { return d & 1 ? portals[d / 2].roomB : portals[d / 2].roomA; };
    for (std::uint32_t d = 0; d < directed && out.listenerRoom != kNone; ++d)
        if (left(d) == out.listenerRoom) portalDistance_[d] = length(portals[d / 2].center - listener_);
    for (std::size_t step = 0; step < directed; ++step) {
        std::uint32_t next = kNone;
        float nearest = INFINITY;
        for (std::uint32_t d = 0; d < directed; ++d) {
            if (!portalDone_[d] && portalDistance_[d] < nearest) {
                nearest = portalDistance_[d];
                next = d;
            }
        }
        if (next == kNone) break;
        portalDone_[next] = 1;
        for (std::uint32_t d = 0; d < directed; ++d) {
            if (portalDone_[d] || left(d) != entered(next)) continue;
            const float through = nearest + length(portals[d / 2].center - portals[next / 2].center);
            if (through < portalDistance_[d]) {
                portalDistance_[d] = through;
                portalPrevious_[d] = next;
            }
        }
    }

    const float blend = 1.f - std::exp(-dt / std::max(config_.smoothingSeconds, 1e-4f));
    out.sourceCount = 0;
    for (std::uint32_t s = 0; s < kMaxSources; ++s) {
        SourceState& source = sources_[s];
        SourceAcoustics& view = out.sources[s];
        view.active = source.active;
        if (!source.active) continue;
        out.sourceCount = s + 1;
        source.room = scene_.roomAt(source.desc.position);
        view.room = source.room;
        view.distance = length(source.desc.position - listener_);
        source.visibility += (source.visibilityTarget - source.visibility) * blend;
        view.visibility = source.visibility;
        for (std::size_t b = 0; b < kBandCount; ++b) {
            source.directEnergy[b] += (source.directTarget[b] - source.directEnergy[b]) * blend;
            view.directGain[b] = std::sqrt(std::max(source.directEnergy[b], 0.f));
        }
        for (std::uint32_t i = 0; i < kMaxPathsPerSource; ++i) {
            PathSlot& slot = source.slots[i];
            ReflectionPath& path = view.reflections[i];
            path = ReflectionPath{};
            if (slot.id == 0) continue;
            for (std::size_t b = 0; b < kBandCount; ++b)
                slot.smoothed[b] += ((slot.visible ? slot.path.amplitude[b] : 0.f) - slot.smoothed[b]) * blend;
            if (!slot.visible && peak(slot.smoothed) < config_.minimumPathGain * 0.1f) {
                slot = PathSlot{};
                continue;
            }
            path.id = slot.id;
            path.delaySeconds = slot.path.length / config_.speedOfSound;
            path.excessSeconds = std::max(slot.path.length - view.distance, 0.f) / config_.speedOfSound;
            path.gain = slot.smoothed;
            path.arrival = slot.path.arrival;
            path.emission = slot.path.emission;
            path.order = slot.key.order;
        }
        updatePortalPaths(s, view);
    }
    published_ = &out;
    exchange_.publish();
}

}  // namespace mspatial::acoustics
