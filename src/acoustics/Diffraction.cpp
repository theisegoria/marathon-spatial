#include "acoustics/Diffraction.hpp"

#include <algorithm>
#include <numbers>
#include <unordered_map>

namespace mspatial::acoustics {
namespace {

constexpr float kPi = std::numbers::pi_v<float>;
constexpr float kWeld = 1e-3f;          // Metres: corners closer than this are one corner.
constexpr float kClearance = 0.02f;     // Validation rays turn this far off the wedge, so they never graze the edge itself.
constexpr float kSectorMargin = 0.17f;  // Radians past half a turn before a sector counts as open air around a convex edge.

Vec3 perpendicular(Vec3 v, Vec3 axis) {
    const Vec3 p = v - axis * dot(v, axis);
    const float n = length(p);
    return n > 1e-9f ? p / n : vec3(0, 0, 0);
}

bool overlaps(const Box& a, const Box& b) {
    return a.min.x <= b.max.x && a.max.x >= b.min.x && a.min.y <= b.max.y && a.max.y >= b.min.y && a.min.z <= b.max.z && a.max.z >= b.min.z;
}

// One edge in world space, with what the search needs of it.
struct WorldEdge {
    EdgeRef ref;
    Vec3 a, axis, h1, h2;
    float length = 0.f;
    bool free = false;  // A free edge or a thin surface: one face, every direction but the face itself is open.
};

// The point on the edge where source -> point -> listener is shortest: where the path unfolds to a straight line about the edge axis.
Vec3 bendPoint(const WorldEdge& e, Vec3 from, Vec3 to) {
    const float sFrom = dot(from - e.a, e.axis), sTo = dot(to - e.a, e.axis);
    const float rFrom = length(from - e.a - e.axis * sFrom), rTo = length(to - e.a - e.axis * sTo);
    const float s = rFrom + rTo > 1e-6f ? (sFrom * rTo + sTo * rFrom) / (rFrom + rTo) : sFrom;
    return e.a + e.axis * std::clamp(s, 0.f, e.length);
}

// Angle of `v` about the edge axis, from h1, in [0, 2 pi).
float angleAbout(const WorldEdge& e, Vec3 v) {
    const Vec3 w = cross(e.axis, e.h1);
    float angle = std::atan2(dot(v, w), dot(v, e.h1));
    if (angle < 0.f) angle += 2.f * kPi;
    return angle;
}

// Which sector between the two faces a point sits in (0 or 1), and how wide that sector is.
std::uint32_t sectorOf(const WorldEdge& e, Vec3 point, Vec3 x, float& width) {
    const float split = angleAbout(e, e.h2);
    const float angle = angleAbout(e, x - point);
    const bool first = angle > 0.f && angle < split;
    width = first ? split : 2.f * kPi - split;
    return first ? 0u : 1u;
}

// True when the point sees the edge as a convex corner: it sits in a sector wider than half a turn.
bool facesCorner(const WorldEdge& e, Vec3 point, Vec3 x) {
    if (e.free) return true;
    float width;
    sectorOf(e, point, x, width);
    return width > kPi + kSectorMargin;
}

// True when the point lies in the plane of one of the faces: a path to it grazes that face, which is fine from either side.
bool onFace(const WorldEdge& e, Vec3 point, Vec3 x) {
    const Vec3 v = x - point - e.axis * dot(x - point, e.axis);
    const float n = length(v);
    if (n < 1e-6f) return true;
    return dot(v / n, e.h1) > 0.9995f || dot(v / n, e.h2) > 0.9995f;
}

// True when both points sit in the same sector, and that sector is open enough for a path to bend around the edge. A point in a face's own
// plane, which a path round the next corner of the same face reaches, counts as being on the other point's side.
bool bendsAround(const WorldEdge& e, Vec3 point, Vec3 x, Vec3 y) {
    if (e.free) return true;
    float wx, wy;
    const std::uint32_t sx = sectorOf(e, point, x, wx), sy = sectorOf(e, point, y, wy);
    if (onFace(e, point, x)) return wy > kPi + kSectorMargin;
    if (onFace(e, point, y)) return wx > kPi + kSectorMargin;
    return sx == sy && wx > kPi + kSectorMargin;
}

Vec3 offWedge(const WorldEdge& e, Vec3 point) {
    const Vec3 inward = e.h1 + e.h2;
    const float n = length(inward);
    return point - (e.free || n < 1e-3f ? e.h1 : inward / n) * kClearance;
}

}  // namespace

void extractDiffractionEdges(Mesh& mesh) {
    mesh.edges.clear();
    const std::size_t triangles = mesh.indices.size() / 3;
    // Weld corners by position, so faces that carry their own copies of a corner still share their edges.
    std::unordered_map<std::uint64_t, std::uint32_t> welded;
    std::vector<std::uint32_t> corner(mesh.vertices.size());
    for (std::uint32_t v = 0; v < mesh.vertices.size(); ++v) {
        const Vec3 p = mesh.vertices[v];
        const auto q = [](float x) { return std::uint64_t(std::int64_t(std::lround(x / kWeld)) & 0x1FFFFF); };
        const std::uint64_t key = q(p.x) | q(p.y) << 21 | q(p.z) << 42;
        corner[v] = welded.try_emplace(key, v).first->second;
    }
    struct Side {
        std::uint32_t opposite;
    };
    std::unordered_map<std::uint64_t, std::vector<Side>> sides;
    for (std::size_t t = 0; t < triangles; ++t) {
        const std::uint32_t idx[3] = {mesh.indices[t * 3], mesh.indices[t * 3 + 1], mesh.indices[t * 3 + 2]};
        const Vec3 a = mesh.vertices[idx[0]], b = mesh.vertices[idx[1]], c = mesh.vertices[idx[2]];
        if (length(cross(b - a, c - a)) < 1e-10f) continue;
        for (std::uint32_t k = 0; k < 3; ++k) {
            const std::uint32_t i = corner[idx[k]], j = corner[idx[(k + 1) % 3]];
            if (i == j) continue;
            const std::uint64_t key = std::uint64_t(std::min(i, j)) << 32 | std::max(i, j);
            sides[key].push_back({idx[(k + 2) % 3]});
        }
    }
    std::vector<Box> boxes;
    for (const auto& [key, list] : sides) {
        const Vec3 a = mesh.vertices[std::uint32_t(key >> 32)], b = mesh.vertices[std::uint32_t(key & 0xFFFFFFFFu)];
        const Vec3 axis = normalize(b - a);
        // Each face's direction away from the edge, within the face.
        std::array<Vec3, 2> directions{};
        std::uint32_t distinct = 0;
        for (const Side& side : list) {
            const Vec3 h = perpendicular(mesh.vertices[side.opposite] - a, axis);
            if (length(h) == 0.f) continue;
            bool known = false;
            for (std::uint32_t d = 0; d < distinct; ++d) known = known || dot(directions[d], h) > 0.999f;
            if (!known && distinct < 2) directions[distinct++] = h;
        }
        if (distinct == 0) continue;
        DiffractionEdge edge;
        edge.a = a;
        edge.b = b;
        edge.h1 = directions[0];
        edge.h2 = distinct == 2 ? directions[1] : directions[0];
        if (distinct == 2 && dot(edge.h1, edge.h2) < -0.999f) continue;  // Flat: one surface carries on across it.
        mesh.edges.push_back(edge);
        Box box;
        box.grow(a);
        box.grow(b);
        box.min = box.min - 0.01f;
        box.max = box.max + 0.01f;
        boxes.push_back(box);
    }
    mesh.edgeBvh.build(boxes);
}

std::uint64_t DiffractionPath::hash() const {
    std::uint64_t h = 0x84222325cbf29ce4ull ^ order;
    for (std::uint32_t i = 0; i < order; ++i) {
        for (std::uint32_t value : {edges[i].instance, edges[i].edge}) {
            h = (h ^ value) * 0x100000001b3ull;
            h ^= h >> 29;
        }
    }
    return h | 1ull;
}

float wedgeLossDb(float detour, float frequency, float speedOfSound) {
    const float fresnel = 2.f * std::max(detour, 0.f) * frequency / speedOfSound;
    const float x = std::sqrt(2.f * kPi * fresnel);
    return 5.f + (x > 1e-4f ? 20.f * std::log10(x / std::tanh(x)) : 0.f);
}

std::uint32_t findDiffractionPaths(const Scene& scene, Vec3 source, Vec3 listener, const DiffractionConfig& config, std::span<DiffractionPath> out) {
    if (out.empty()) return 0;
    const float direct = length(listener - source);
    if (direct < 1e-3f) return 0;
    // Every point of a path no longer than direct + maxDetour lies in an ellipsoid; this box holds it.
    const float grow = 0.5f * std::sqrt(2.f * direct * config.maxDetour + config.maxDetour * config.maxDetour) + 0.01f;
    Box query;
    query.grow(source);
    query.grow(listener);
    query.min = query.min - grow;
    query.max = query.max + grow;

    // Convex corners facing the source or the listener, shortest detour first. A path over one edge needs both to see the same corner; a
    // path over two needs the source to see the first and the listener the second, so either is worth keeping for pairing.
    struct Candidate {
        WorldEdge edge;
        float detour;
        bool single;
    };
    std::vector<Candidate> singles;
    const auto& instances = scene.instances();
    for (std::uint32_t id = 0; id < instances.size(); ++id) {
        const Instance& instance = instances[id];
        if (!instance.live || !instance.enabled || !overlaps(instance.worldBounds, query)) continue;
        const Mesh& mesh = scene.meshes()[instance.mesh];
        if (mesh.edges.empty() || mesh.edgeBvh.nodes.empty()) continue;
        Box local;
        for (int k = 0; k < 8; ++k)
            local.grow(instance.worldToObject.point(vec3(k & 1 ? query.max.x : query.min.x, k & 2 ? query.max.y : query.min.y, k & 4 ? query.max.z : query.min.z)));
        std::uint32_t stack[64];
        int top = 0;
        stack[top++] = 0;
        while (top > 0) {
            const BvhNode& node = mesh.edgeBvh.nodes[stack[--top]];
            if (!overlaps(node.bounds, local)) continue;
            if (node.count == 0) {
                if (top + 2 > 64) continue;
                stack[top++] = node.first;
                stack[top++] = node.first + 1;
                continue;
            }
            for (std::uint32_t i = 0; i < node.count; ++i) {
                const std::uint32_t index = mesh.edgeBvh.order[node.first + i];
                const DiffractionEdge& objectEdge = mesh.edges[index];
                WorldEdge e;
                e.ref = {id, index};
                e.a = instance.objectToWorld.point(objectEdge.a);
                const Vec3 b = instance.objectToWorld.point(objectEdge.b);
                e.length = length(b - e.a);
                if (e.length < 1e-4f) continue;
                e.axis = (b - e.a) / e.length;
                e.h1 = perpendicular(instance.objectToWorld.direction(objectEdge.h1), e.axis);
                e.h2 = perpendicular(instance.objectToWorld.direction(objectEdge.h2), e.axis);
                e.free = dot(e.h1, e.h2) > 0.999f;
                const Vec3 p = bendPoint(e, source, listener);
                const float detour = length(p - source) + length(listener - p) - direct;
                if (detour < 1e-4f || detour > config.maxDetour) continue;
                if (!facesCorner(e, p, source) && !facesCorner(e, p, listener)) continue;
                singles.push_back({e, detour, bendsAround(e, p, source, listener)});
            }
        }
    }
    std::sort(singles.begin(), singles.end(), [](const Candidate& x, const Candidate& y) { return x.detour < y.detour; });
    if (singles.size() > config.candidates) singles.resize(config.candidates);

    const auto finish = [&](DiffractionPath& path, const std::array<float, kMaxDiffractionOrder>& detours) {
        const Vec3 last = path.points[path.order - 1];
        path.arrival = normalize(last - listener);
        path.emission = normalize(path.points[0] - source);
        for (std::uint32_t b = 0; b < kBandCount; ++b) {
            float loss = 0.f;
            for (std::uint32_t k = 0; k < path.order; ++k) loss += wedgeLossDb(detours[k], kBandCenterHz[b], config.speedOfSound);
            path.amplitude[b] = std::pow(10.f, -loss / 20.f) * std::exp(-0.5f * config.air[b] * path.length) / std::max(path.length, 1.f);
        }
    };
    std::vector<DiffractionPath> found;
    for (const Candidate& c : singles) {
        if (!c.single) continue;
        DiffractionPath path;
        path.order = 1;
        path.edges[0] = c.edge.ref;
        path.points[0] = bendPoint(c.edge, source, listener);
        path.clear[0] = offWedge(c.edge, path.points[0]);
        path.length = direct + c.detour;
        finish(path, {c.detour, 0.f});
        found.push_back(path);
    }
    // Two edges, around both corners of something thick. Each bend point is placed by turns until the pair settles.
    for (std::size_t i = 0; i < singles.size(); ++i) {
        for (std::size_t j = 0; j < singles.size(); ++j) {
            if (i == j) continue;
            const WorldEdge& e1 = singles[i].edge;
            const WorldEdge& e2 = singles[j].edge;
            Vec3 p1 = bendPoint(e1, source, listener), p2 = bendPoint(e2, source, listener);
            for (int pass = 0; pass < 4; ++pass) {
                p1 = bendPoint(e1, source, p2);
                p2 = bendPoint(e2, p1, listener);
            }
            const float a = length(p1 - source), m = length(p2 - p1), z = length(listener - p2);
            if (m < 0.05f) continue;  // The same corner twice.
            const float total = a + m + z;
            if (total - direct > config.maxDetour) continue;
            if (!bendsAround(e1, p1, source, p2) || !bendsAround(e2, p2, p1, listener)) continue;
            const float d1 = a + m - length(p2 - source), d2 = m + z - length(listener - p1);
            if (d1 < 1e-4f || d2 < 1e-4f) continue;
            DiffractionPath path;
            path.order = 2;
            path.edges = {e1.ref, e2.ref};
            path.points = {p1, p2};
            path.clear = {offWedge(e1, p1), offWedge(e2, p2)};
            path.length = total;
            finish(path, {d1, d2});
            found.push_back(path);
        }
    }
    std::sort(found.begin(), found.end(), [](const DiffractionPath& x, const DiffractionPath& y) { return x.amplitude[4] > y.amplitude[4]; });
    const std::uint32_t count = std::min<std::uint32_t>(std::uint32_t(found.size()), std::uint32_t(out.size()));
    std::copy_n(found.begin(), count, out.begin());
    return count;
}

}  // namespace mspatial::acoustics
