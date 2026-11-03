#include "acoustics/AcousticScene.hpp"
#include "acoustics/Diffraction.hpp"

#include <algorithm>
#include <cassert>
#include <unordered_map>

namespace mspatial::acoustics {
namespace {

constexpr int kBins = 12;
constexpr std::uint32_t kLeafSize = 4;

Vec3 centroid(const Box& b) { return (b.min + b.max) * 0.5f; }

std::uint32_t buildNode(Bvh& bvh, std::span<const Box> boxes, std::uint32_t nodeIndex, std::uint32_t begin, std::uint32_t end) {
    Box bounds, centers;
    for (std::uint32_t i = begin; i < end; ++i) {
        bounds.grow(boxes[bvh.order[i]]);
        centers.grow(centroid(boxes[bvh.order[i]]));
    }
    bvh.nodes[nodeIndex].bounds = bounds;
    const std::uint32_t count = end - begin;
    const auto leaf = [&] {
        bvh.nodes[nodeIndex].first = begin;
        bvh.nodes[nodeIndex].count = count;
        return nodeIndex;
    };
    if (count <= kLeafSize) return leaf();

    float bestCost = INFINITY;
    int bestAxis = -1;
    float bestSplit = 0.f;
    const Vec3 extent = centers.max - centers.min;
    for (int axis = 0; axis < 3; ++axis) {
        if (extent[axis] <= 1e-7f) continue;
        Box binBox[kBins];
        std::uint32_t binCount[kBins] = {};
        const float scale = float(kBins) / extent[axis];
        const auto binOf = [&](const Box& b) {
            return std::min(kBins - 1, int((centroid(b)[axis] - centers.min[axis]) * scale));
        };
        for (std::uint32_t i = begin; i < end; ++i) {
            const Box& b = boxes[bvh.order[i]];
            const int bin = binOf(b);
            binBox[bin].grow(b);
            ++binCount[bin];
        }
        float rightArea[kBins] = {};
        std::uint32_t rightCount[kBins] = {};
        Box accum;
        std::uint32_t total = 0;
        for (int bin = kBins - 1; bin > 0; --bin) {
            accum.grow(binBox[bin]);
            total += binCount[bin];
            rightArea[bin] = accum.area();
            rightCount[bin] = total;
        }
        accum = Box{};
        total = 0;
        for (int bin = 0; bin < kBins - 1; ++bin) {
            accum.grow(binBox[bin]);
            total += binCount[bin];
            if (total == 0 || rightCount[bin + 1] == 0) continue;
            const float cost = accum.area() * float(total) + rightArea[bin + 1] * float(rightCount[bin + 1]);
            if (cost < bestCost) {
                bestCost = cost;
                bestAxis = axis;
                bestSplit = centers.min[axis] + float(bin + 1) / scale;
            }
        }
    }
    if (bestAxis < 0) return leaf();
    const auto middle = std::partition(bvh.order.begin() + begin, bvh.order.begin() + end,
                                       [&](std::uint32_t i) { return centroid(boxes[i])[bestAxis] < bestSplit; });
    const std::uint32_t mid = std::uint32_t(middle - bvh.order.begin());
    if (mid == begin || mid == end) return leaf();
    const std::uint32_t left = std::uint32_t(bvh.nodes.size());
    bvh.nodes.emplace_back();
    bvh.nodes.emplace_back();
    bvh.nodes[nodeIndex].first = left;
    bvh.nodes[nodeIndex].count = 0;
    buildNode(bvh, boxes, left, begin, mid);
    buildNode(bvh, boxes, left + 1, mid, end);
    return nodeIndex;
}

Transform inverseAffine(const Transform& t) {
    const float (*m)[4] = t.rows;
    const float c00 = m[1][1] * m[2][2] - m[1][2] * m[2][1], c01 = m[0][2] * m[2][1] - m[0][1] * m[2][2],
                c02 = m[0][1] * m[1][2] - m[0][2] * m[1][1];
    const float c10 = m[1][2] * m[2][0] - m[1][0] * m[2][2], c11 = m[0][0] * m[2][2] - m[0][2] * m[2][0],
                c12 = m[0][2] * m[1][0] - m[0][0] * m[1][2];
    const float c20 = m[1][0] * m[2][1] - m[1][1] * m[2][0], c21 = m[0][1] * m[2][0] - m[0][0] * m[2][1],
                c22 = m[0][0] * m[1][1] - m[0][1] * m[1][0];
    const float det = m[0][0] * c00 + m[0][1] * c10 + m[0][2] * c20;
    const float inv = 1.f / det;
    Transform out;
    const float linear[3][3] = {{c00 * inv, c01 * inv, c02 * inv}, {c10 * inv, c11 * inv, c12 * inv}, {c20 * inv, c21 * inv, c22 * inv}};
    for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 3; ++c) out.rows[r][c] = linear[r][c];
        out.rows[r][3] = -(linear[r][0] * m[0][3] + linear[r][1] * m[1][3] + linear[r][2] * m[2][3]);
    }
    return out;
}

Transform normalMatrix(const Transform& inverse) {
    Transform out;
    for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 3; ++c) out.rows[r][c] = inverse.rows[c][r];
        out.rows[r][3] = 0.f;
    }
    return out;
}

void updateWorld(Instance& instance, const Mesh& mesh) {
    instance.worldToObject = inverseAffine(instance.objectToWorld);
    instance.normalToWorld = normalMatrix(instance.worldToObject);
    instance.worldBounds = Box{};
    for (int corner = 0; corner < 8; ++corner) {
        const Vec3 p = vec3((corner & 1) ? mesh.bounds.max.x : mesh.bounds.min.x, (corner & 2) ? mesh.bounds.max.y : mesh.bounds.min.y,
                            (corner & 4) ? mesh.bounds.max.z : mesh.bounds.min.z);
        instance.worldBounds.grow(instance.objectToWorld.point(p));
    }
}

bool slab(const Box& box, Vec3 origin, Vec3 inverse, float tMin, float tMax) {
    const Vec3 a = (box.min - origin) * inverse;
    const Vec3 b = (box.max - origin) * inverse;
    const Vec3 near = simd_min(a, b);
    const Vec3 far = simd_max(a, b);
    const float enter = std::max({near.x, near.y, near.z, tMin});
    const float leave = std::min({far.x, far.y, far.z, tMax});
    return enter <= leave;
}

Vec3 safeInverse(Vec3 d) {
    const auto one = [](float v) { return 1.f / (std::fabs(v) < 1e-12f ? std::copysign(1e-12f, v) : v); };
    return vec3(one(d.x), one(d.y), one(d.z));
}

bool triangle(Vec3 origin, Vec3 direction, Vec3 a, Vec3 b, Vec3 c, float tMin, float& tBest) {
    const Vec3 e1 = b - a, e2 = c - a;
    const Vec3 p = cross(direction, e2);
    const float det = dot(e1, p);
    if (std::fabs(det) < 1e-12f) return false;
    const float inv = 1.f / det;
    const Vec3 s = origin - a;
    const float u = dot(s, p) * inv;
    if (u < 0.f || u > 1.f) return false;
    const Vec3 q = cross(s, e1);
    const float v = dot(direction, q) * inv;
    if (v < 0.f || u + v > 1.f) return false;
    const float t = dot(e2, q) * inv;
    if (t < tMin || t >= tBest) return false;
    tBest = t;
    return true;
}

template <class Leaf>
void traverse(const Bvh& bvh, Vec3 origin, Vec3 inverse, float tMin, float& tMax, Leaf&& leaf) {
    if (bvh.nodes.empty()) return;
    std::uint32_t stack[64];
    int top = 0;
    stack[top++] = 0;
    while (top > 0) {
        const BvhNode& node = bvh.nodes[stack[--top]];
        if (!slab(node.bounds, origin, inverse, tMin, tMax)) continue;
        if (node.count != 0) {
            for (std::uint32_t i = 0; i < node.count; ++i) leaf(bvh.order[node.first + i]);
        } else {
            stack[top++] = node.first;
            stack[top++] = node.first + 1;
        }
    }
}

}  // namespace

void Bvh::build(std::span<const Box> boxes) {
    nodes.clear();
    order.resize(boxes.size());
    for (std::uint32_t i = 0; i < boxes.size(); ++i) order[i] = i;
    if (boxes.empty()) return;
    nodes.reserve(boxes.size() * 2);
    nodes.emplace_back();
    buildNode(*this, boxes, 0, 0, std::uint32_t(boxes.size()));
}

std::uint32_t Scene::addMaterial(const Material& material) {
    materials_.push_back(material);
    ++materialRevision_;
    return std::uint32_t(materials_.size() - 1);
}

void Scene::setMaterial(std::uint32_t id, const Material& material) {
    materials_.at(id) = material;
    ++materialRevision_;
}

std::uint32_t Scene::addMesh(const MeshDesc& desc) {
    assert(desc.positions.size() % 3 == 0 && desc.indices.size() % 3 == 0);
    assert(desc.triangleMaterials.size() == desc.indices.size() / 3);
    Mesh mesh;
    mesh.live = true;
    mesh.vertices.reserve(desc.positions.size() / 3);
    for (std::size_t i = 0; i + 2 < desc.positions.size(); i += 3) {
        mesh.vertices.push_back(vec3(desc.positions[i], desc.positions[i + 1], desc.positions[i + 2]));
        mesh.bounds.grow(mesh.vertices.back());
    }
    mesh.indices.assign(desc.indices.begin(), desc.indices.end());
    const std::size_t triangles = mesh.indices.size() / 3;
    mesh.triangleFacet.assign(triangles, 0);

    std::unordered_map<std::uint64_t, std::uint32_t> planes;
    std::vector<Box> boxes(triangles);
    for (std::size_t t = 0; t < triangles; ++t) {
        const Vec3 a = mesh.vertices[mesh.indices[t * 3]], b = mesh.vertices[mesh.indices[t * 3 + 1]],
                   c = mesh.vertices[mesh.indices[t * 3 + 2]];
        boxes[t].grow(a);
        boxes[t].grow(b);
        boxes[t].grow(c);
        const Vec3 raw = cross(b - a, c - a);
        const float area = length(raw);
        if (area < 1e-10f) continue;  // Degenerate: rays never hit it, and it stays on facet 0.
        Vec3 normal = raw / area;
        // A plane has two faces and one reflection, so the sign of the winding is dropped.
        const bool flip = normal.x < -1e-6f || (std::fabs(normal.x) <= 1e-6f && (normal.y < -1e-6f || (std::fabs(normal.y) <= 1e-6f && normal.z < 0.f)));
        if (flip) normal = -normal;
        const float offset = dot(normal, a);
        const auto quantize = [](float v, float scale) { return std::uint64_t(std::int64_t(std::lround(v * scale)) & 0x7FFF); };
        const std::uint64_t key = quantize(normal.x, 2048.f) | quantize(normal.y, 2048.f) << 15 | quantize(normal.z, 2048.f) << 30 |
                                  (std::uint64_t(std::int64_t(std::lround(offset * 256.f))) & 0x3FFF) << 45 |
                                  std::uint64_t(desc.triangleMaterials[t] & 31) << 59;
        const auto [entry, inserted] = planes.try_emplace(key, std::uint32_t(mesh.facets.size()));
        if (inserted) {
            Facet facet;
            facet.normal = normal;
            facet.offset = offset;
            facet.material = desc.triangleMaterials[t];
            mesh.facets.push_back(facet);
        }
        mesh.triangleFacet[t] = entry->second;
        mesh.facets[entry->second].bounds.grow(boxes[t]);
    }
    if (mesh.facets.empty()) mesh.facets.emplace_back();
    mesh.bvh.build(boxes);
    extractDiffractionEdges(mesh);

    ++meshRevision_;
    mesh.serial = meshRevision_;
    if (!freeMeshes_.empty()) {
        const std::uint32_t id = freeMeshes_.back();
        freeMeshes_.pop_back();
        meshes_[id] = std::move(mesh);
        return id;
    }
    meshes_.push_back(std::move(mesh));
    return std::uint32_t(meshes_.size() - 1);
}

void Scene::removeMesh(std::uint32_t id) {
    for (std::uint32_t i = 0; i < instances_.size(); ++i)
        if (instances_[i].live && instances_[i].mesh == id) removeInstance(i);
    meshes_.at(id) = Mesh{};
    freeMeshes_.push_back(id);
    ++meshRevision_;
}

std::uint32_t Scene::addInstance(std::uint32_t mesh, const Transform& transform, bool enabled) {
    Instance instance;
    instance.mesh = mesh;
    instance.objectToWorld = transform;
    instance.enabled = enabled;
    instance.live = true;
    updateWorld(instance, meshes_.at(mesh));
    ++instanceRevision_;
    if (!freeInstances_.empty()) {
        const std::uint32_t id = freeInstances_.back();
        freeInstances_.pop_back();
        instances_[id] = instance;
        return id;
    }
    instances_.push_back(instance);
    return std::uint32_t(instances_.size() - 1);
}

void Scene::setInstanceTransform(std::uint32_t id, const Transform& transform) {
    Instance& instance = instances_.at(id);
    instance.objectToWorld = transform;
    updateWorld(instance, meshes_[instance.mesh]);
    ++instanceRevision_;
}

void Scene::setInstanceEnabled(std::uint32_t id, bool enabled) {
    instances_.at(id).enabled = enabled;
    ++instanceRevision_;
}

void Scene::removeInstance(std::uint32_t id) {
    instances_.at(id) = Instance{};
    freeInstances_.push_back(id);
    ++instanceRevision_;
}

std::uint32_t Scene::addRoom(const RoomDesc& desc) {
    rooms_.push_back(desc);
    return std::uint32_t(rooms_.size() - 1);
}

std::uint32_t Scene::addPortal(const PortalDesc& desc) {
    portals_.push_back(desc);
    return std::uint32_t(portals_.size() - 1);
}

void Scene::setPortal(std::uint32_t id, const PortalDesc& desc) { portals_.at(id) = desc; }

bool Scene::facetValid(std::uint32_t instance, std::uint32_t facet) const {
    if (instance >= instances_.size()) return false;
    const Instance& i = instances_[instance];
    return i.live && i.enabled && facet < meshes_[i.mesh].facets.size();
}

const Facet& Scene::facet(std::uint32_t instance, std::uint32_t facet) const {
    return meshes_[instances_[instance].mesh].facets[facet];
}

const Material& Scene::facetMaterial(std::uint32_t instance, std::uint32_t id) const { return materials_[facet(instance, id).material]; }

WorldPlane Scene::worldPlane(std::uint32_t instance, std::uint32_t id) const {
    const Instance& i = instances_[instance];
    const Facet& f = facet(instance, id);
    WorldPlane out;
    out.normal = normalize(i.normalToWorld.direction(f.normal));
    out.offset = dot(out.normal, i.objectToWorld.point(f.normal * f.offset));
    return out;
}

bool Scene::nearFacet(std::uint32_t instance, std::uint32_t id, Vec3 world, float margin) const {
    return facet(instance, id).bounds.contains(instances_[instance].worldToObject.point(world), margin);
}

std::uint32_t Scene::roomAt(Vec3 point) const {
    std::uint32_t best = kNone;
    float bestVolume = INFINITY;
    for (std::uint32_t r = 0; r < rooms_.size(); ++r) {
        float volume = 0.f;
        bool inside = false;
        for (const Box& box : rooms_[r].bounds) {
            volume += box.volume();
            inside = inside || box.contains(point);
        }
        if (inside && volume < bestVolume) {
            best = r;
            bestVolume = volume;
        }
    }
    return best;
}

void SceneIntersector::rebuild(const Scene& scene) {
    enabled_.clear();
    std::vector<Box> boxes;
    const auto& instances = scene.instances();
    for (std::uint32_t i = 0; i < instances.size(); ++i) {
        if (!instances[i].live || !instances[i].enabled) continue;
        enabled_.push_back(i);
        boxes.push_back(instances[i].worldBounds);
    }
    top_.build(boxes);
}

RayHit SceneIntersector::closest(const Scene& scene, Vec3 origin, Vec3 direction, float tMin, float tMax) const {
    RayHit best;
    float tBest = tMax;
    traverse(top_, origin, safeInverse(direction), tMin, tBest, [&](std::uint32_t slot) {
        const std::uint32_t id = enabled_[slot];
        const Instance& instance = scene.instances()[id];
        const Mesh& mesh = scene.meshes()[instance.mesh];
        // The object-space direction keeps its scale, so t stays in world units.
        const Vec3 objectOrigin = instance.worldToObject.point(origin);
        const Vec3 objectDirection = instance.worldToObject.direction(direction);
        traverse(mesh.bvh, objectOrigin, safeInverse(objectDirection), tMin, tBest, [&](std::uint32_t t) {
            const Vec3 a = mesh.vertices[mesh.indices[t * 3]], b = mesh.vertices[mesh.indices[t * 3 + 1]],
                       c = mesh.vertices[mesh.indices[t * 3 + 2]];
            if (triangle(objectOrigin, objectDirection, a, b, c, tMin, tBest)) {
                best.instance = id;
                best.facet = mesh.triangleFacet[t];
                best.t = tBest;
            }
        });
    });
    return best;
}

}  // namespace mspatial::acoustics
