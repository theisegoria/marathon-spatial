#include "acoustics/EnergyField.hpp"

#include <algorithm>
#include <numbers>

namespace mspatial::acoustics {
namespace {

constexpr float kPi = std::numbers::pi_v<float>;
constexpr float kMaxSegment = 200.f;  // Metres per ray: past this a path is in the open and ends.

std::array<Vec3, kFieldDirections> makeDirections() {
    const float t = (1.f + std::sqrt(5.f)) / 2.f;
    const float raw[kFieldDirections][3] = {{-1, t, 0}, {1, t, 0}, {-1, -t, 0}, {1, -t, 0}, {0, -1, t}, {0, 1, t},
                                            {0, -1, -t}, {0, 1, -t}, {t, 0, -1}, {t, 0, 1}, {-t, 0, -1}, {-t, 0, 1}};
    std::array<Vec3, kFieldDirections> out;
    for (std::uint32_t i = 0; i < kFieldDirections; ++i) out[i] = normalize(vec3(raw[i][0], raw[i][1], raw[i][2]));
    return out;
}

float uniform(std::uint32_t& state) {
    state ^= state << 13;
    state ^= state >> 17;
    state ^= state << 5;
    return float(state >> 8) * (1.f / 16777216.f);
}

// Cosine-weighted direction about a unit normal.
Vec3 lambert(Vec3 normal, float u, float v) {
    const float r = std::sqrt(u), phi = 2.f * kPi * v;
    const Vec3 helper = std::fabs(normal.z) < 0.9f ? vec3(0, 0, 1) : vec3(1, 0, 0);
    const Vec3 t = normalize(cross(normal, helper)), b = cross(normal, t);
    return normalize(t * (r * std::cos(phi)) + b * (r * std::sin(phi)) + normal * std::sqrt(std::max(0.f, 1.f - u)));
}

}  // namespace

const std::array<Vec3, kFieldDirections>& fieldDirections() {
    static const std::array<Vec3, kFieldDirections> directions = makeDirections();
    return directions;
}

std::uint32_t nearestFieldDirection(Vec3 direction) {
    const auto& directions = fieldDirections();
    std::uint32_t best = 0;
    float bestDot = -2.f;
    for (std::uint32_t i = 0; i < kFieldDirections; ++i) {
        const float d = dot(directions[i], direction);
        if (d > bestDot) {
            bestDot = d;
            best = i;
        }
    }
    return best;
}

float EnergyField::band(std::uint32_t b, std::uint32_t bin) const {
    float sum = 0.f;
    for (std::uint32_t d = 0; d < kFieldDirections; ++d) sum += at(d, b, bin);
    return sum;
}

Bands EnergyField::t60() const {
    Bands out = uniformBands(0.f);
    if (bins < 4) return out;
    std::vector<double> edc(bins);
    for (std::uint32_t b = 0; b < kBandCount; ++b) {
        double running = 0.0;
        for (std::uint32_t i = bins; i-- > 0;) {
            running += band(b, i);
            edc[i] = running;
        }
        if (edc[0] <= 0.0) continue;
        // Least squares on the decay curve between -5 dB and the deepest point the field supports, capped at -35 dB. The last tenth is left
        // out: a truncated integral drops toward zero there and would steepen the fit.
        double sx = 0, sy = 0, sxx = 0, sxy = 0;
        std::uint32_t n = 0;
        float deepest = 0.f;
        const std::uint32_t usable = bins - bins / 10;
        for (std::uint32_t i = 0; i < usable; ++i) {
            if (edc[i] <= 0.0) break;
            const float db = 10.f * std::log10(float(edc[i] / edc[0]));
            if (db > -5.f) continue;
            if (db < -35.f) break;
            const double x = (double(i) + 0.5) * binSeconds;
            sx += x;
            sy += db;
            sxx += x * x;
            sxy += x * db;
            ++n;
            deepest = db;
        }
        if (n < 3 || deepest > -15.f) continue;
        const double slope = (double(n) * sxy - sx * sy) / (double(n) * sxx - sx * sx);
        if (slope < 0.0) out[b] = float(-60.0 / slope);
    }
    return out;
}

void EnergyTracer::configure(const EnergyTracerConfig& config) {
    config_ = config;
    bins_ = std::max<std::uint32_t>(1, std::uint32_t(std::ceil(config.seconds / config.binSeconds)));
    radius_ = config.receiverRadius;
    paths_.assign(config.paths, Path{});
    sums_.assign(std::size_t(kFieldDirections) * kBandCount * bins_, 0.0);
    emitted_ = 0.0;
    ++generation_;
}

void EnergyTracer::setRoomVolume(float volume) {
    // A sphere of about 0.15 V^(1/3) catches the same fraction of a diffuse field in any room of the same shape.
    radius_ = volume > 0.f ? std::max(config_.receiverRadius, 0.15f * std::cbrt(volume)) : config_.receiverRadius;
}

void EnergyTracer::setRoom(std::uint32_t room) {
    if (room == room_) return;
    room_ = room;
    reset();
}

void EnergyTracer::reset() {
    std::fill(sums_.begin(), sums_.end(), 0.0);
    emitted_ = 0.0;
    ++generation_;
    for (Path& path : paths_) {
        path.live = false;
        path.inFlight = false;
    }
}

void EnergyTracer::age(float dt) {
    if (config_.timeConstant <= 0.f || dt <= 0.f) return;
    const double keep = std::exp(-double(dt) / config_.timeConstant);
    for (double& value : sums_) value *= keep;
    emitted_ *= keep;
}

void EnergyTracer::spawn(Path& path, Vec3 listener) {
    // One R2 sequence across all spawns, so the sphere of initial directions fills evenly however paths retire.
    const std::uint32_t index = sequence_++;
    const float u = float(index * 3242174889u) / 4294967296.f, v = float(index * 2447445414u) / 4294967296.f;
    const float z = 1.f - 2.f * u, ring = std::sqrt(std::max(0.f, 1.f - z * z)), phi = 2.f * kPi * v;
    path.origin = listener;
    path.direction = vec3(ring * std::cos(phi), ring * std::sin(phi), z);
    path.travelled = 0.f;
    path.energy = uniformBands(1.f);
    path.order = 0;
    path.generation = generation_;
    path.live = true;
    emitted_ += 1.0;
}

std::uint32_t EnergyTracer::issue(TraceBatch& batch, std::uint32_t budget, Vec3 listener, float selfHitOffset, std::vector<EnergyJob>& jobs) {
    std::uint32_t issued = 0;
    const float horizon = config_.speedOfSound * config_.seconds;
    for (std::uint32_t p = 0; p < paths_.size() && issued < budget; ++p) {
        Path& path = paths_[p];
        if (path.inFlight) continue;
        if (!path.live || path.generation != generation_) spawn(path, listener);
        const float reach = std::min(kMaxSegment, horizon - path.travelled);
        if (reach <= 0.f) {
            path.live = false;
            continue;
        }
        jobs.push_back({p, batch.closestCount, path.generation});
        batch.closestRays[batch.closestCount++] = {path.origin.x, path.origin.y, path.origin.z, path.order == 0 ? 0.f : selfHitOffset,
                                                   path.direction.x, path.direction.y, path.direction.z, reach};
        path.inFlight = true;
        ++issued;
    }
    return issued;
}

void EnergyTracer::consume(const Scene& scene, const TraceBatch& batch, std::span<const EnergyJob> jobs, Vec3 listener) {
    const float volume = 4.f / 3.f * kPi * radius_ * radius_ * radius_;
    const float scale = 4.f * kPi / volume;
    for (const EnergyJob& job : jobs) {
        if (job.path >= paths_.size()) continue;
        Path& path = paths_[job.path];
        path.inFlight = false;
        if (job.generation != path.generation || !path.live) continue;
        if (batch.failed) continue;  // The segment is simply asked again.
        ++segments_;
        const AcousticHitRecord& hit = batch.closestHits[job.ray];
        const AcousticRayRecord& ray = batch.closestRays[job.ray];
        const bool missed = hit.instance == ACOUSTIC_MISS || !scene.facetValid(hit.instance, hit.facet);
        const float length = missed ? ray.tMax : hit.t;

        // Score the segment against the receiver sphere. The segment leaving the listener is the direct sound and does not count.
        if (path.order > 0) {
            const Vec3 offset = path.origin - listener;
            const float b = dot(path.direction, offset), c = dot(offset, offset) - radius_ * radius_;
            const float disc = b * b - c;
            if (disc > 0.f) {
                const float root = std::sqrt(disc);
                const float s0 = std::max(0.f, -b - root), s1 = std::min(length, -b + root);
                if (s1 > s0) {
                    const float middle = 0.5f * (s0 + s1);
                    const float arrival = (path.travelled + middle) / config_.speedOfSound;
                    const auto bin = std::uint32_t(arrival / config_.binSeconds);
                    if (bin < bins_) {
                        // The sound travels along the path, so it arrives from the opposite way.
                        const std::uint32_t d = nearestFieldDirection(-path.direction);
                        const float chord = (s1 - s0) * scale;
                        for (std::uint32_t band = 0; band < kBandCount; ++band) {
                            const float energy = path.energy[band] * std::exp(-config_.air[band] * middle);
                            sums_[(std::size_t(d) * kBandCount + band) * bins_ + bin] += double(energy * chord);
                        }
                    }
                }
            }
        }

        if (missed) {
            path.live = false;  // Into the open: nothing comes back.
            continue;
        }
        path.travelled += length;
        const Material& material = scene.facetMaterial(hit.instance, hit.facet);
        float strongest = 0.f, scattering = 0.f;
        for (std::uint32_t band = 0; band < kBandCount; ++band) {
            path.energy[band] *= (1.f - std::clamp(material.absorption[band], 0.f, 1.f)) * std::exp(-config_.air[band] * length);
            strongest = std::max(strongest, path.energy[band]);
            scattering += material.scattering[band] / float(kBandCount);
        }
        if (path.travelled >= config_.speedOfSound * config_.seconds || strongest < config_.energyFloor) {
            path.live = false;
            continue;
        }
        Vec3 normal = scene.worldPlane(hit.instance, hit.facet).normal;
        if (dot(normal, path.direction) > 0.f) normal = -normal;
        path.origin = path.origin + path.direction * length;
        if (uniform(rng_) < scattering)
            path.direction = lambert(normal, uniform(rng_), uniform(rng_));
        else
            path.direction = normalize(path.direction - normal * (2.f * dot(path.direction, normal)));
        ++path.order;
    }
}

void EnergyTracer::read(EnergyField& out) const {
    out.bins = bins_;
    out.binSeconds = config_.binSeconds;
    out.room = room_;
    out.paths = emitted_;
    out.energy.resize(sums_.size());
    const double norm = emitted_ > 0.0 ? 1.0 / emitted_ : 0.0;
    for (std::size_t i = 0; i < sums_.size(); ++i) out.energy[i] = float(sums_[i] * norm);
}

}  // namespace mspatial::acoustics
