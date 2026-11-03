// The live acoustic scene, owned by one non-real-time thread. Meshes are static shapes; instances place them and move every tick.
// Triangles are grouped into facets, one per plane and material, because an image source belongs to a plane: a wall made of many
// triangles is one reflector, and a reflection walking across its diagonal keeps its identity.
#pragma once

#include "acoustics/AcousticTypes.hpp"

namespace mspatial::acoustics {

struct BvhNode {
    Box bounds;
    std::uint32_t first = 0;  // First child, or first primitive of a leaf.
    std::uint32_t count = 0;  // Primitive count; 0 for an interior node.
};

// Binned surface-area heuristic over boxes. `order` maps leaf primitive ranges to input indices.
struct Bvh {
    std::vector<BvhNode> nodes;
    std::vector<std::uint32_t> order;
    void build(std::span<const Box> boxes);
};

struct Facet {
    Vec3 normal = vec3(0, 0, 1);  // Object space, unit length.
    float offset = 0.f;           // dot(normal, x) == offset on the plane.
    std::uint32_t material = 0;
    Box bounds;
};

// An edge that can diffract: where two faces meet at an angle, or a free edge of a thin surface. Object space. h1 and h2 are unit vectors
// perpendicular to the edge, lying in each face and pointing away from the edge into it; a free edge has h2 == h1.
struct DiffractionEdge {
    Vec3 a = vec3(0, 0, 0), b = vec3(0, 0, 0);
    Vec3 h1 = vec3(1, 0, 0), h2 = vec3(1, 0, 0);
};

struct Mesh {
    std::vector<Vec3> vertices;
    std::vector<std::uint32_t> indices;
    std::vector<std::uint32_t> triangleFacet;
    std::vector<Facet> facets;
    Bvh bvh;
    Box bounds;
    std::uint64_t serial = 0;  // Unique per addMesh, so a reused slot reads as a new mesh.
    std::vector<DiffractionEdge> edges;
    Bvh edgeBvh;
    bool live = false;
};

struct Instance {
    std::uint32_t mesh = kNone;
    Transform objectToWorld;
    Transform worldToObject;
    Transform normalToWorld;  // Inverse transpose of the linear part.
    Box worldBounds;
    bool enabled = true;
    bool live = false;
};

struct WorldPlane {
    Vec3 normal = vec3(0, 0, 1);
    float offset = 0.f;
};

class Scene {
public:
    std::uint32_t addMaterial(const Material& material);
    void setMaterial(std::uint32_t id, const Material& material);

    std::uint32_t addMesh(const MeshDesc& desc);
    void removeMesh(std::uint32_t id);  // Removes instances of it too.

    std::uint32_t addInstance(std::uint32_t mesh, const Transform& transform, bool enabled = true);
    void setInstanceTransform(std::uint32_t id, const Transform& transform);
    void setInstanceEnabled(std::uint32_t id, bool enabled);
    void removeInstance(std::uint32_t id);

    std::uint32_t addRoom(const RoomDesc& desc);
    std::uint32_t addPortal(const PortalDesc& desc);
    void setPortal(std::uint32_t id, const PortalDesc& desc);

    const std::vector<Material>& materials() const { return materials_; }
    const std::vector<Mesh>& meshes() const { return meshes_; }
    const std::vector<Instance>& instances() const { return instances_; }
    const std::vector<RoomDesc>& rooms() const { return rooms_; }
    const std::vector<PortalDesc>& portals() const { return portals_; }

    // Changes when triangles appear or vanish.
    std::uint64_t meshRevision() const { return meshRevision_; }
    // Changes when anything rays can see moves, including instance transforms and enabled flags.
    std::uint64_t instanceRevision() const { return instanceRevision_; }
    std::uint64_t materialRevision() const { return materialRevision_; }

    bool facetValid(std::uint32_t instance, std::uint32_t facet) const;
    const Facet& facet(std::uint32_t instance, std::uint32_t facet) const;
    const Material& facetMaterial(std::uint32_t instance, std::uint32_t facet) const;
    WorldPlane worldPlane(std::uint32_t instance, std::uint32_t facet) const;
    // True when the world point is inside the facet's object-space bounds grown by margin.
    bool nearFacet(std::uint32_t instance, std::uint32_t facet, Vec3 world, float margin) const;

    // The smallest room whose bounds contain the point, or kNone.
    std::uint32_t roomAt(Vec3 point) const;

private:
    std::vector<Material> materials_;
    std::vector<Mesh> meshes_;
    std::vector<Instance> instances_;
    std::vector<RoomDesc> rooms_;
    std::vector<PortalDesc> portals_;
    std::vector<std::uint32_t> freeMeshes_;
    std::vector<std::uint32_t> freeInstances_;
    std::uint64_t meshRevision_ = 1;
    std::uint64_t instanceRevision_ = 1;
    std::uint64_t materialRevision_ = 1;
};

// Closest hit of a world-space ray over a top-level BVH of the enabled instances. Shared by the CPU backend and tests.
struct RayHit {
    std::uint32_t instance = kNone;
    std::uint32_t facet = 0;
    float t = 0.f;
};

class SceneIntersector {
public:
    void rebuild(const Scene& scene);
    RayHit closest(const Scene& scene, Vec3 origin, Vec3 direction, float tMin, float tMax) const;

private:
    Bvh top_;
    std::vector<std::uint32_t> enabled_;
};

}  // namespace mspatial::acoustics
