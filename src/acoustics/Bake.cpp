#include "acoustics/Bake.hpp"

#include "acoustics/AcousticTrace.hpp"

#include <algorithm>
#include <atomic>
#include <bit>
#include <cstring>
#include <fstream>
#include <numbers>
#include <thread>
#include <unordered_set>

namespace mspatial::acoustics {
namespace {

static_assert(std::endian::native == std::endian::little, "bake files are written as the machine holds them");

constexpr float kPi = std::numbers::pi_v<float>;
constexpr std::uint32_t kNearest = 4;
constexpr std::uint32_t kWindows = 2;
constexpr float kShadowRamp = 20.f * kPi / 180.f;
constexpr float kLineEps = 1e-3f;
constexpr char kMagic[8] = {'M', 'S', 'B', 'A', 'K', 'E', '\r', '\n'};
constexpr std::uint32_t kVersion = 1;
constexpr std::uint32_t kVolumeRays = 128;
constexpr float kVolumeRange = 40.f;
constexpr std::uint32_t kMaxRoute = 64;

struct Fnv {
    std::uint64_t h = 0xcbf29ce484222325ull;
    void bytes(const void* data, std::size_t size) {
        const auto* c = static_cast<const unsigned char*>(data);
        for (std::size_t i = 0; i < size; ++i) h = (h ^ c[i]) * 0x100000001b3ull;
    }
    void value(float v) { bytes(&v, sizeof v); }
    void value(std::uint32_t v) { bytes(&v, sizeof v); }
    void value(Vec3 v) {
        value(v.x);
        value(v.y);
        value(v.z);
    }
    void value(const Bands& b) {
        for (float v : b) value(v);
    }
};

// Runs body(state, i) for every i below count on up to `threads` threads, the calling one among them. Each thread makes its own state with
// make. onCaller runs on the calling thread after each item it finishes, with how many are done in all.
template <class Make, class Body>
void parallelFor(std::uint32_t count, std::uint32_t threads, const Make& make, const Body& body, const std::function<void(std::uint32_t)>& onCaller) {
    if (count == 0) return;
    std::atomic<std::uint32_t> cursor{0}, done{0};
    const auto work = [&](bool caller) {
        auto state = make();
        for (;;) {
            const std::uint32_t i = cursor.fetch_add(1, std::memory_order_relaxed);
            if (i >= count) break;
            body(state, i);
            const std::uint32_t finished = done.fetch_add(1, std::memory_order_acq_rel) + 1;
            if (caller && onCaller) onCaller(finished);
        }
    };
    std::vector<std::thread> pool;
    const std::uint32_t workers = std::min(std::max(threads, 1u), count);
    for (std::uint32_t t = 1; t < workers; ++t) pool.emplace_back(work, false);
    work(true);
    for (std::thread& t : pool) t.join();
}

Vec3 fibonacci(std::uint32_t index, std::uint32_t count) {
    const float z = 1.f - (2.f * float(index) + 1.f) / float(count);
    const float r = std::sqrt(std::max(0.f, 1.f - z * z));
    const float phi = 2.f * kPi * std::fmod(float(index) * 0.61803398875f, 1.f);
    return vec3(r * std::cos(phi), r * std::sin(phi), z);
}

// The 26 directions to the neighbours of a cube cell: faces, edges and corners.
std::array<Vec3, 26> neighbourDirections() {
    std::array<Vec3, 26> out;
    std::uint32_t n = 0;
    for (int x = -1; x <= 1; ++x)
        for (int y = -1; y <= 1; ++y)
            for (int z = -1; z <= 1; ++z)
                if (x || y || z) out[n++] = normalize(vec3(float(x), float(y), float(z)));
    return out;
}

// Two parallel rays 2 mm apart must both pass. A single ray can slip between two triangles exactly along their shared edge, and on a grid of
// probes some line always does: one such leak would join the probes inside a solid to the room outside.
bool clearLine(const Scene& scene, const SceneIntersector& rays, Vec3 a, Vec3 b) {
    const Vec3 d = b - a;
    const float span = length(d);
    if (span < 2.f * kLineEps) return true;
    const Vec3 dir = d / span;
    if (rays.closest(scene, a, dir, kLineEps, span - kLineEps).instance != kNone) return false;
    const Vec3 helper = std::fabs(dir.y) < 0.9f ? vec3(0, 1, 0) : vec3(1, 0, 0);
    const Vec3 side = normalize(cross(dir, helper) + helper * 0.37f) * 2e-3f;
    return rays.closest(scene, a + side, dir, kLineEps, span - kLineEps).instance == kNone;
}

// V = 4 pi / 3 <r^3> over the sphere of directions, the estimate the simulator's room probes make.
float estimateVolume(const Scene& scene, const SceneIntersector& rays, Vec3 point) {
    double sum = 0.0;
    for (std::uint32_t i = 0; i < kVolumeRays; ++i) {
        const RayHit hit = rays.closest(scene, point, fibonacci(i, kVolumeRays), 0.f, kVolumeRange);
        const double r = hit.instance == kNone ? kVolumeRange : hit.t;
        sum += r * r * r;
    }
    return float(4.0 / 3.0 * std::numbers::pi * sum / kVolumeRays);
}

// Energy of input bins (width inSeconds) falling in output bin j (width outSeconds).
float rebinned(const float* in, std::uint32_t inBins, float inSeconds, std::uint32_t j, float outSeconds) {
    const double t0 = double(j) * outSeconds, t1 = t0 + outSeconds;
    const auto k0 = std::uint32_t(std::max(0.0, std::floor(t0 / inSeconds)));
    const auto k1 = std::uint32_t(std::min<double>(inBins, std::ceil(t1 / inSeconds)));
    double sum = 0.0;
    for (std::uint32_t k = k0; k < k1; ++k) {
        const double a = std::max(t0, double(k) * inSeconds), b = std::min(t1, double(k + 1) * inSeconds);
        if (b > a) sum += double(in[k]) * (b - a) / inSeconds;
    }
    return float(sum);
}

template <class T>
void put(std::ofstream& file, const T& value) {
    file.write(reinterpret_cast<const char*>(&value), sizeof value);
}
template <class T>
void putAll(std::ofstream& file, const std::vector<T>& values) {
    file.write(reinterpret_cast<const char*>(values.data()), std::streamsize(values.size() * sizeof(T)));
}
template <class T>
bool get(std::ifstream& file, T& value) {
    return bool(file.read(reinterpret_cast<char*>(&value), sizeof value));
}
template <class T>
bool getAll(std::ifstream& file, std::vector<T>& values, std::size_t count) {
    values.resize(count);
    return bool(file.read(reinterpret_cast<char*>(values.data()), std::streamsize(count * sizeof(T))));
}

struct Phases {
    const std::function<void(float)>* progress;
    float start = 0.f, span = 0.f;
    void begin(float share) {
        start += span;
        span = share;
        report(0.f);
    }
    void report(float done) const {
        if (*progress) (*progress)(std::clamp(start + span * done, 0.f, 1.f));
    }
};

}  // namespace

bool AcousticBake::sees(std::uint32_t a, std::uint32_t b) const {
    const std::uint32_t words = (probeCount() + 63) / 64;
    return (visibility[std::size_t(a) * words + b / 64] >> (b % 64)) & 1u;
}

bool AcousticBake::fieldAt(Vec3 point, std::uint32_t room, std::uint32_t outBins, float outSeconds, EnergyField& out) const {
    if (!hasReverb() || outBins == 0 || !(outSeconds > 0.f)) return false;
    std::array<std::uint32_t, kNearest> ids{};
    std::array<float, kNearest> d2{};
    std::uint32_t found = 0;
    for (std::uint32_t p = 0; p < probeCount(); ++p) {
        if (rooms[p] != room) continue;
        const Vec3 offset = positions[p] - point;
        const float d = dot(offset, offset);
        std::uint32_t at = found < kNearest ? found++ : kNearest;
        if (at == kNearest) {
            if (d >= d2[kNearest - 1]) continue;
            at = kNearest - 1;
        }
        while (at > 0 && d2[at - 1] > d) {
            ids[at] = ids[at - 1];
            d2[at] = d2[at - 1];
            --at;
        }
        ids[at] = p;
        d2[at] = d;
    }
    if (found == 0) return false;
    std::array<float, kNearest> weight{};
    float sum = 0.f;
    for (std::uint32_t k = 0; k < found; ++k) sum += weight[k] = 1.f / (d2[k] + 0.25f);
    out.bins = outBins;
    out.binSeconds = outSeconds;
    out.room = room;
    out.paths = pathsPerProbe;
    out.energy.assign(std::size_t(kFieldDirections) * kBandCount * outBins, 0.f);
    for (std::uint32_t k = 0; k < found; ++k) {
        const std::uint32_t p = ids[k];
        const float w = weight[k] / sum;
        const float* share = shares.data() + std::size_t(p) * kWindows * kFieldDirections * kBandCount;
        for (std::uint32_t b = 0; b < kBandCount; ++b) {
            const float* total = totals.data() + (std::size_t(p) * kBandCount + b) * bins;
            for (std::uint32_t j = 0; j < outBins; ++j) {
                const float e = rebinned(total, bins, binSeconds, j, outSeconds) * w;
                if (e == 0.f) continue;
                const std::uint32_t window = float(j) * outSeconds < earlySeconds ? 0 : 1;
                for (std::uint32_t d = 0; d < kFieldDirections; ++d) out.at(d, b, j) += e * share[(window * kFieldDirections + d) * kBandCount + b];
            }
        }
    }
    return true;
}

void finishProbePath(DiffractionPath& path, Vec3 source, Vec3 listener, const DiffractionConfig& config) {
    std::array<Vec3, kMaxBendPoints + 2> p;
    p[0] = source;
    for (std::uint32_t k = 0; k < path.order; ++k) p[k + 1] = path.points[k];
    p[path.order + 1] = listener;
    path.length = 0.f;
    for (std::uint32_t k = 0; k <= path.order; ++k) path.length += length(p[k + 1] - p[k]);
    const auto unit = [](Vec3 v) {
        const float n = length(v);
        return n > 1e-6f ? v / n : vec3(0, 0, 0);
    };
    path.emission = unit(p[1] - p[0]);
    path.arrival = unit(p[path.order] - listener);
    // A probe stands in open air, away from the corner it stands in for, so a detour measured at each probe overstates every bend. The route
    // as a whole is a better measure: one wedge loss for its whole detour over the straight line, with the 5 dB at the shadow boundary
    // ramped in as the sharpest bend reaches 20 degrees.
    float sharpest = 0.f;
    for (std::uint32_t k = 0; k < path.order; ++k) {
        const float turn = std::acos(std::clamp(dot(unit(p[k + 1] - p[k]), unit(p[k + 2] - p[k + 1])), -1.f, 1.f));
        sharpest = std::max(sharpest, turn);
    }
    const float detour = std::max(0.f, path.length - length(listener - source));
    const float boundary = 5.f * std::min(sharpest / kShadowRamp, 1.f);
    for (std::uint32_t b = 0; b < kBandCount; ++b) {
        const float loss = wedgeLossDb(detour, kBandCenterHz[b], config.speedOfSound) - 5.f + boundary;
        path.amplitude[b] = std::pow(10.f, -loss / 20.f) * std::exp(-0.5f * config.air[b] * path.length) / std::max(path.length, 1.f);
    }
}

std::uint32_t AcousticBake::findPaths(Vec3 source, Vec3 listener, const DiffractionConfig& config, std::uint32_t rotate,
                                      std::span<DiffractionPath> out) const {
    const std::uint32_t n = probeCount();
    if (!hasPathing() || out.empty() || n == 0) return 0;
    const float direct = length(listener - source);
    const auto nearest = [&](Vec3 point, std::array<std::uint32_t, kNearest>& ids) {
        std::array<float, kNearest> d2{};
        std::uint32_t found = 0;
        for (std::uint32_t p = 0; p < n; ++p) {
            const Vec3 offset = positions[p] - point;
            const float d = dot(offset, offset);
            std::uint32_t at = found < kNearest ? found++ : kNearest;
            if (at == kNearest) {
                if (d >= d2[kNearest - 1]) continue;
                at = kNearest - 1;
            }
            while (at > 0 && d2[at - 1] > d) {
                ids[at] = ids[at - 1];
                d2[at] = d2[at - 1];
                --at;
            }
            ids[at] = p;
            d2[at] = d;
        }
        return found;
    };
    std::array<std::uint32_t, kNearest> fromSource{}, toListener{};
    const std::uint32_t ns = nearest(source, fromSource), nl = nearest(listener, toListener);
    std::vector<DiffractionPath> found;
    std::array<std::uint32_t, kMaxRoute> route{}, tight{};
    for (std::uint32_t i = 0; i < ns; ++i) {
        for (std::uint32_t j = 0; j < nl; ++j) {
            const std::uint32_t a = fromSource[i], b = toListener[j];
            if (a != b && !std::isfinite(distances[std::size_t(a) * n + b])) continue;
            std::uint32_t count = 0, at = a;
            route[count++] = a;
            while (at != b && count < kMaxRoute) route[count++] = at = next[std::size_t(at) * n + b];
            if (at != b) continue;
            // Pull the route tight: from each probe, on to the furthest one along it that it sees.
            std::uint32_t m = 0;
            tight[m++] = route[0];
            for (std::uint32_t k = 0; k + 1 < count;) {
                std::uint32_t far = count - 1;
                while (far > k + 1 && !sees(route[k], route[far])) --far;
                tight[m++] = route[far];
                k = far;
            }
            // The source may well see past the first probe and the listener past the last, so routes without them are candidates too.
            for (std::uint32_t skipFirst : {1u, 0u}) {
                for (std::uint32_t skipLast : {1u, 0u}) {
                    if (m <= skipFirst + skipLast) continue;
                    const std::uint32_t bends = m - skipFirst - skipLast;
                    if (bends > kMaxBendPoints) continue;
                    DiffractionPath path;
                    path.order = bends;
                    for (std::uint32_t k = 0; k < bends; ++k) {
                        const std::uint32_t probe = tight[skipFirst + k];
                        path.edges[k] = {kBakeProbeTag, probe};
                        path.points[k] = path.clear[k] = positions[probe];
                    }
                    finishProbePath(path, source, listener, config);
                    if (path.length - direct > config.maxRouteDetour) continue;
                    const std::uint64_t id = path.hash();
                    if (std::none_of(found.begin(), found.end(), [&](const DiffractionPath& f) { return f.hash() == id; })) found.push_back(path);
                }
            }
        }
    }
    if (found.empty()) return 0;
    std::sort(found.begin(), found.end(), [](const DiffractionPath& x, const DiffractionPath& y) { return x.amplitude[4] > y.amplitude[4]; });
    const auto count = std::uint32_t(std::min(out.size(), found.size()));
    const std::size_t base = std::size_t(rotate) * count;
    for (std::uint32_t k = 0; k < count; ++k) out[k] = found[(base + k) % found.size()];
    return count;
}

std::uint64_t staticSceneHash(const Scene& scene) {
    Fnv fnv;
    const auto& instances = scene.instances();
    for (std::uint32_t id = 0; id < instances.size(); ++id) {
        const Instance& instance = instances[id];
        if (!instance.live || !instance.enabled) continue;
        const Mesh& mesh = scene.meshes()[instance.mesh];
        fnv.value(id);
        fnv.bytes(instance.objectToWorld.rows, sizeof instance.objectToWorld.rows);
        for (Vec3 v : mesh.vertices) fnv.value(v);
        fnv.bytes(mesh.indices.data(), mesh.indices.size() * sizeof(std::uint32_t));
        for (std::uint32_t facet : mesh.triangleFacet) {
            const Material& m = scene.materials()[mesh.facets[facet].material];
            fnv.value(m.absorption);
            fnv.value(m.scattering);
            fnv.value(m.transmission);
        }
    }
    return fnv.h;
}

bool bakeScene(const Scene& scene, const BakeConfig& config, AcousticBake& out, std::string* error, const std::function<void(float)>& progress) {
    out = AcousticBake{};
    const auto fail = [&](const char* what) {
        if (error) *error = what;
        out = AcousticBake{};
        return false;
    };
    const bool grid = config.probes.empty();
    if (grid && !(config.spacing >= 0.1f)) return fail("probe spacing must be at least 10 cm");
    if (config.reverb && (!(config.energy.seconds > 0.f) || !(config.energy.binSeconds > 0.f) || config.pathsPerProbe == 0))
        return fail("the reverb bake needs a field length, a bin width and paths");
    const std::uint32_t threads = std::clamp(config.threads, 1u, 64u);
    Phases phases{&progress};
    SceneIntersector rays;
    rays.rebuild(scene);

    // Probes.
    phases.begin(0.05f);
    std::vector<Vec3> placed;
    std::vector<std::uint32_t> placedRooms;
    std::uint32_t dropped = 0;
    const auto& sceneRooms = scene.rooms();
    if (!grid) {
        if (config.probes.size() > kMaxBakeProbes) return fail("more than 4096 probes");
        placed = config.probes;
        for (Vec3 p : placed) placedRooms.push_back(scene.roomAt(p));
    } else {
        struct Region {
            Box box;
            std::uint32_t room;
        };
        std::vector<Region> regions;
        for (std::uint32_t r = 0; r < sceneRooms.size(); ++r)
            for (const Box& box : sceneRooms[r].bounds) regions.push_back({box, r});
        if (regions.empty()) {
            Box all;
            for (const Instance& instance : scene.instances())
                if (instance.live && instance.enabled) all.grow(instance.worldBounds);
            if (all.min.x <= all.max.x) regions.push_back({all, kNone});
        }
        std::unordered_set<std::uint64_t> taken;
        const float cell = 0.5f * config.spacing;
        const auto directions = neighbourDirections();
        for (const Region& region : regions) {
            const Vec3 extent = simd_max(region.box.max - region.box.min, vec3(0, 0, 0));
            const auto steps = [&](float e) { return std::max<std::uint64_t>(1, std::uint64_t(std::lround(e / config.spacing))); };
            const std::uint64_t nx = steps(extent.x), ny = steps(extent.y), nz = steps(extent.z);
            if (nx * ny * nz > 16ull * kMaxBakeProbes) return fail("the probe grid would hold far more than 4096 probes: raise the spacing");
            for (std::uint64_t ix = 0; ix < nx; ++ix)
                for (std::uint64_t iy = 0; iy < ny; ++iy)
                    for (std::uint64_t iz = 0; iz < nz; ++iz) {
                        const Vec3 p = region.box.min + extent * vec3((float(ix) + 0.5f) / float(nx), (float(iy) + 0.5f) / float(ny),
                                                                      (float(iz) + 0.5f) / float(nz));
                        // Where rooms overlap, the smallest owns the point.
                        const std::uint32_t room = scene.roomAt(p);
                        if (!sceneRooms.empty() && room != region.room) continue;
                        const auto quantize = [&](float v) { return std::uint64_t(std::int64_t(std::floor(v / cell)) & 0x1fffff); };
                        if (!taken.insert(quantize(p.x) | quantize(p.y) << 21 | quantize(p.z) << 42).second) continue;
                        bool clear = true;
                        for (Vec3 d : directions)
                            if (rays.closest(scene, p, d, 0.f, config.clearance).instance != kNone) {
                                clear = false;
                                break;
                            }
                        if (!clear) {
                            ++dropped;
                            continue;
                        }
                        placed.push_back(p);
                        placedRooms.push_back(room);
                    }
        }
    }
    if (placed.empty()) return fail("no probe stands clear of the geometry");

    // Which probes see each other: every grid bake needs it to drop enclosed probes, and pathing is built on it.
    const auto n0 = std::uint32_t(placed.size());
    std::vector<std::uint8_t> seen;
    const bool needVisibility = config.pathing || grid;
    phases.begin(needVisibility ? 0.15f : 0.f);
    if (needVisibility) {
        if (n0 > 2 * kMaxBakeProbes) return fail("more than 8192 candidate probes: raise the spacing");
        seen.assign(std::size_t(n0) * n0, 0);
        parallelFor(
            n0, threads,
            [&] {
                SceneIntersector local;
                local.rebuild(scene);
                return local;
            },
            [&](SceneIntersector& local, std::uint32_t a) {
                for (std::uint32_t b = a + 1; b < n0; ++b) seen[std::size_t(a) * n0 + b] = clearLine(scene, local, placed[a], placed[b]);
            },
            [&](std::uint32_t done) { phases.report(float(done) / float(n0)); });
        for (std::uint32_t a = 0; a < n0; ++a)
            for (std::uint32_t b = a + 1; b < n0; ++b) seen[std::size_t(b) * n0 + a] = seen[std::size_t(a) * n0 + b];
    }

    // A grid probe in a group that sees no other group, and is much smaller than the largest, stands inside something solid.
    std::vector<std::uint32_t> index(n0, kNone);
    {
        std::vector<std::uint32_t> parent(n0);
        for (std::uint32_t i = 0; i < n0; ++i) parent[i] = i;
        const auto root = [&](std::uint32_t i) {
            while (parent[i] != i) i = parent[i] = parent[parent[i]];
            return i;
        };
        std::vector<std::uint32_t> size(n0, 0);
        std::uint32_t minimum = 0;
        if (grid) {
            for (std::uint32_t a = 0; a < n0; ++a)
                for (std::uint32_t b = a + 1; b < n0; ++b)
                    if (seen[std::size_t(a) * n0 + b]) parent[root(a)] = root(b);
            std::uint32_t largest = 0;
            for (std::uint32_t i = 0; i < n0; ++i) largest = std::max(largest, ++size[root(i)]);
            minimum = std::min(std::max(2u, (largest + 9) / 10), largest);
        }
        std::uint32_t kept = 0;
        for (std::uint32_t i = 0; i < n0; ++i) {
            if (grid && size[root(i)] < minimum) {
                ++dropped;
                continue;
            }
            index[i] = kept++;
            out.positions.push_back(placed[i]);
            out.rooms.push_back(placedRooms[i]);
        }
    }
    const std::uint32_t n = out.probeCount();
    if (n > kMaxBakeProbes) return fail("more than 4096 probes: raise the spacing");
    out.sceneHash = staticSceneHash(scene);
    out.stats.probes = n;
    out.stats.droppedProbes = dropped;

    // Routes: the shortest way between every pair through probes that see each other.
    phases.begin(config.pathing ? 0.10f : 0.f);
    if (config.pathing) {
        const std::uint32_t words = (n + 63) / 64;
        out.visibility.assign(std::size_t(n) * words, 0);
        for (std::uint32_t a = 0; a < n0; ++a)
            for (std::uint32_t b = 0; b < n0; ++b)
                if (a != b && index[a] != kNone && index[b] != kNone && seen[std::size_t(a) * n0 + b]) {
                    out.visibility[std::size_t(index[a]) * words + index[b] / 64] |= 1ull << (index[b] % 64);
                    out.stats.visiblePairs += index[a] < index[b];
                }
        out.distances.assign(std::size_t(n) * n, INFINITY);
        out.next.assign(std::size_t(n) * n, 0);
        struct Scratch {
            std::vector<float> distance;
            std::vector<std::uint8_t> settled;
            std::vector<std::uint16_t> first;
        };
        parallelFor(
            n, threads, [&] { return Scratch{std::vector<float>(n), std::vector<std::uint8_t>(n), std::vector<std::uint16_t>(n)}; },
            [&](Scratch& s, std::uint32_t a) {
                std::fill(s.distance.begin(), s.distance.end(), INFINITY);
                std::fill(s.settled.begin(), s.settled.end(), 0);
                s.distance[a] = 0.f;
                s.first[a] = std::uint16_t(a);
                for (;;) {
                    std::uint32_t u = kNone;
                    float best = INFINITY;
                    for (std::uint32_t v = 0; v < n; ++v)
                        if (!s.settled[v] && s.distance[v] < best) {
                            best = s.distance[v];
                            u = v;
                        }
                    if (u == kNone) break;
                    s.settled[u] = 1;
                    const std::uint64_t* row = out.visibility.data() + std::size_t(u) * words;
                    for (std::uint32_t w = 0; w < words; ++w) {
                        for (std::uint64_t bits = row[w]; bits; bits &= bits - 1) {
                            const std::uint32_t v = w * 64 + std::uint32_t(std::countr_zero(bits));
                            if (s.settled[v]) continue;
                            const float through = best + length(out.positions[v] - out.positions[u]);
                            if (through < s.distance[v]) {
                                s.distance[v] = through;
                                s.first[v] = u == a ? std::uint16_t(v) : s.first[u];
                            }
                        }
                    }
                }
                std::copy(s.distance.begin(), s.distance.end(), out.distances.begin() + std::ptrdiff_t(a) * n);
                std::copy(s.first.begin(), s.first.end(), out.next.begin() + std::ptrdiff_t(a) * n);
                for (std::uint32_t v = 0; v < n; ++v)
                    if (!std::isfinite(s.distance[v])) out.next[std::size_t(a) * n + v] = std::uint16_t(a);
            },
            [&](std::uint32_t done) { phases.report(float(done) / float(n)); });
    }

    // Reverb: each probe's field, traced to the end with every path counted only once it has finished.
    phases.begin(config.reverb ? 0.70f : 0.f);
    if (config.reverb) {
        EnergyTracerConfig energy = config.energy;
        energy.paths = std::max(energy.paths, 1u);
        energy.timeConstant = 0.f;
        EnergyTracer probe;
        probe.configure(energy);
        EnergyField shape;
        probe.read(shape);
        out.bins = shape.bins;
        out.binSeconds = energy.binSeconds;
        out.pathsPerProbe = config.pathsPerProbe;
        out.totals.assign(std::size_t(n) * kBandCount * out.bins, 0.f);
        out.shares.assign(std::size_t(n) * kWindows * kFieldDirections * kBandCount, 0.f);
        std::atomic<std::uint64_t> segments{0};
        struct Tracing {
            std::unique_ptr<TraceBackend> backend;
            std::unique_ptr<EnergyTracer> tracer;
            SceneIntersector rays;
            std::vector<EnergyJob> jobs;
            EnergyField field;
        };
        parallelFor(
            n, threads,
            [&] {
                Tracing t;
                TraceCapacity capacity;
                capacity.closestRays = energy.paths;
                capacity.layeredRays = 1;
                capacity.discoveryRays = 1;
                t.backend = makeCpuTraceBackend(capacity);
                t.backend->syncScene(scene);
                t.tracer = std::make_unique<EnergyTracer>();
                t.tracer->configure(energy);
                t.rays.rebuild(scene);
                return t;
            },
            [&](Tracing& t, std::uint32_t p) {
                EnergyTracer& tracer = *t.tracer;
                const Vec3 at = out.positions[p];
                const std::uint32_t room = out.rooms[p];
                tracer.reset();
                tracer.setSpawning(true);
                tracer.setRoomVolume(room != kNone && sceneRooms[room].volume > 0.f ? sceneRooms[room].volume : estimateVolume(scene, t.rays, at));
                const std::uint64_t before = tracer.segmentsTraced();
                bool spawning = true;
                for (std::uint32_t guard = 0; guard < (1u << 24); ++guard) {
                    if (spawning && tracer.emitted() >= double(config.pathsPerProbe)) tracer.setSpawning(spawning = false);
                    if (!spawning && tracer.livePaths() == 0) break;
                    TraceBatch* batch = t.backend->acquireBatch();
                    if (!batch) break;
                    t.jobs.clear();
                    tracer.issue(*batch, energy.paths, at, t.backend->selfHitOffset(), t.jobs);
                    t.backend->submit(batch);
                    TraceBatch* done = t.backend->pollCompleted();
                    if (done) {
                        tracer.consume(scene, *done, t.jobs, at);
                        t.backend->release(done);
                    }
                }
                segments += tracer.segmentsTraced() - before;
                tracer.read(t.field);
                const EnergyField& field = t.field;
                for (std::uint32_t b = 0; b < kBandCount; ++b) {
                    float* total = out.totals.data() + (std::size_t(p) * kBandCount + b) * out.bins;
                    std::array<std::array<double, kFieldDirections>, kWindows> windowed{};
                    for (std::uint32_t i = 0; i < field.bins; ++i) {
                        const std::uint32_t window = float(i) * field.binSeconds < out.earlySeconds ? 0 : 1;
                        double sum = 0.0;
                        for (std::uint32_t d = 0; d < kFieldDirections; ++d) {
                            sum += field.at(d, b, i);
                            windowed[window][d] += field.at(d, b, i);
                        }
                        total[i] = float(sum);
                    }
                    for (std::uint32_t w = 0; w < kWindows; ++w) {
                        double all = 0.0;
                        for (double v : windowed[w]) all += v;
                        float* share = out.shares.data() + (std::size_t(p) * kWindows + w) * kFieldDirections * kBandCount;
                        for (std::uint32_t d = 0; d < kFieldDirections; ++d)
                            share[d * kBandCount + b] = all > 0.0 ? float(windowed[w][d] / all) : 1.f / float(kFieldDirections);
                    }
                }
            },
            [&](std::uint32_t done) { phases.report(float(done) / float(n)); });
        out.stats.segments = segments.load();
    }
    if (progress) progress(1.f);
    return true;
}

void blendFields(const EnergyField& live, const EnergyField& baked, double priorPaths, EnergyField& out) {
    const double n = std::max(live.paths, 0.0), prior = std::max(priorPaths, 0.0), total = n + prior;
    out.bins = live.bins;
    out.binSeconds = live.binSeconds;
    out.room = live.room;
    out.paths = total;
    out.energy.resize(live.energy.size());
    const bool fits = baked.bins == live.bins && baked.energy.size() == live.energy.size();
    const float wl = total > 0.0 ? float(n / total) : 1.f, wb = fits && total > 0.0 ? float(prior / total) : 0.f;
    for (std::size_t i = 0; i < live.energy.size(); ++i) out.energy[i] = wl * live.energy[i] + (fits ? wb * baked.energy[i] : 0.f);
    if (!fits) out.paths = n;
}

bool AcousticBake::save(const std::string& path, std::string* error) const {
    std::ofstream file(path, std::ios::binary | std::ios::trunc);
    if (!file) {
        if (error) *error = "cannot open the file for writing";
        return false;
    }
    const std::uint32_t n = probeCount(), flags = (hasReverb() ? 1u : 0u) | (hasPathing() ? 2u : 0u);
    file.write(kMagic, sizeof kMagic);
    put(file, kVersion);
    put(file, n);
    put(file, flags);
    put(file, sceneHash);
    put(file, bins);
    put(file, binSeconds);
    put(file, earlySeconds);
    put(file, pathsPerProbe);
    put(file, stats.probes);
    put(file, stats.droppedProbes);
    put(file, stats.visiblePairs);
    put(file, stats.segments);
    for (Vec3 p : positions) {
        put(file, p.x);
        put(file, p.y);
        put(file, p.z);
    }
    putAll(file, rooms);
    if (hasReverb()) {
        putAll(file, totals);
        putAll(file, shares);
    }
    if (hasPathing()) {
        putAll(file, visibility);
        putAll(file, distances);
        putAll(file, next);
    }
    if (!file.flush()) {
        if (error) *error = "could not write the whole file";
        return false;
    }
    return true;
}

bool AcousticBake::load(const std::string& path, AcousticBake& out, std::string* error) {
    out = AcousticBake{};
    const auto fail = [&](const char* what) {
        if (error) *error = what;
        out = AcousticBake{};
        return false;
    };
    std::ifstream file(path, std::ios::binary);
    if (!file) return fail("cannot open the file");
    char magic[sizeof kMagic];
    std::uint32_t version = 0, n = 0, flags = 0;
    if (!file.read(magic, sizeof magic) || std::memcmp(magic, kMagic, sizeof magic) != 0) return fail("not a Marathon Spatial bake");
    if (!get(file, version) || version != kVersion) return fail("a bake from another version");
    if (!get(file, n) || !get(file, flags) || !get(file, out.sceneHash) || !get(file, out.bins) || !get(file, out.binSeconds) ||
        !get(file, out.earlySeconds) || !get(file, out.pathsPerProbe) || !get(file, out.stats.probes) || !get(file, out.stats.droppedProbes) ||
        !get(file, out.stats.visiblePairs) || !get(file, out.stats.segments))
        return fail("the header is cut short");
    if (n == 0 || n > kMaxBakeProbes || (flags & ~3u) != 0) return fail("the header is out of range");
    const bool reverb = flags & 1u, pathing = flags & 2u;
    if (reverb && (out.bins == 0 || out.bins > 100000 || !(out.binSeconds > 0.f) || !(out.pathsPerProbe > 0.0)))
        return fail("the reverb header is out of range");
    if (!reverb) out.bins = 0;
    out.positions.resize(n);
    for (Vec3& p : out.positions) {
        float x, y, z;
        if (!get(file, x) || !get(file, y) || !get(file, z)) return fail("the probes are cut short");
        p = vec3(x, y, z);
    }
    if (!getAll(file, out.rooms, n)) return fail("the probes are cut short");
    if (reverb && (!getAll(file, out.totals, std::size_t(n) * kBandCount * out.bins) ||
                   !getAll(file, out.shares, std::size_t(n) * kWindows * kFieldDirections * kBandCount)))
        return fail("the reverb data is cut short");
    if (pathing) {
        const std::size_t words = (n + 63) / 64;
        if (!getAll(file, out.visibility, n * words) || !getAll(file, out.distances, std::size_t(n) * n) || !getAll(file, out.next, std::size_t(n) * n))
            return fail("the routes are cut short");
        for (std::uint16_t v : out.next)
            if (v >= n) return fail("a route leads to a probe that does not exist");
        for (float d : out.distances)
            if (!(d >= 0.f)) return fail("a route has a negative length");
    }
    if (file.peek() != std::ifstream::traits_type::eof()) return fail("the file is longer than its contents");
    return true;
}

}  // namespace mspatial::acoustics
