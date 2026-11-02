// Marathon Spatial acoustic ray queries. Three kernels over one instance acceleration structure: closest hits for image-source
// validation, room probes and portal openness; up to ACOUSTIC_MAX_LAYERS hits per ray for transmission; and listener discovery walks.
#include <metal_raytracing>
#include "AcousticTraceShared.h"

using namespace metal;
using namespace metal::raytracing;

typedef intersector<triangle_data, instancing> AcousticIntersector;

static AcousticIntersector makeIntersector() {
    AcousticIntersector query;
    query.assume_geometry_type(geometry_type::triangle);
    query.force_opacity(forced_opacity::opaque);
    query.accept_any_intersection(false);
    query.set_triangle_cull_mode(triangle_cull_mode::none);
    return query;
}

static AcousticHitRecord record(AcousticIntersector::result_type hit, const device AcousticInstanceRecord* instances,
                                const device uint* triangleFacets) {
    AcousticHitRecord out;
    out.instance = ACOUSTIC_MISS;
    out.facet = 0;
    out.t = hit.distance;
    out.pad = 0;
    if (hit.type == intersection_type::none) return out;
    out.instance = hit.user_instance_id;
    out.facet = triangleFacets[instances[out.instance].firstTriangle + hit.primitive_id];
    return out;
}

kernel void acoustic_closest(constant AcousticTraceParams& params [[buffer(0)]],
                             const device AcousticRayRecord* rays [[buffer(1)]],
                             device AcousticHitRecord* hits [[buffer(2)]],
                             const device AcousticInstanceRecord* instances [[buffer(3)]],
                             const device uint* triangleFacets [[buffer(4)]],
                             instance_acceleration_structure scene [[buffer(5)]],
                             uint index [[thread_position_in_grid]]) {
    const AcousticRayRecord source = rays[index];
    ray query;
    query.origin = float3(source.originX, source.originY, source.originZ);
    query.direction = float3(source.directionX, source.directionY, source.directionZ);
    query.min_distance = source.tMin;
    query.max_distance = source.tMax;
    AcousticIntersector intersect = makeIntersector();
    hits[index] = record(intersect.intersect(query, scene), instances, triangleFacets);
}

kernel void acoustic_layered(constant AcousticTraceParams& params [[buffer(0)]],
                             const device AcousticRayRecord* rays [[buffer(1)]],
                             device AcousticLayeredHit* hits [[buffer(2)]],
                             const device AcousticInstanceRecord* instances [[buffer(3)]],
                             const device uint* triangleFacets [[buffer(4)]],
                             instance_acceleration_structure scene [[buffer(5)]],
                             uint index [[thread_position_in_grid]]) {
    const AcousticRayRecord source = rays[index];
    AcousticLayeredHit out;
    out.count = 0;
    out.pad0 = out.pad1 = out.pad2 = 0;
    ray query;
    query.origin = float3(source.originX, source.originY, source.originZ);
    query.direction = float3(source.directionX, source.directionY, source.directionZ);
    query.min_distance = source.tMin;
    query.max_distance = source.tMax;
    AcousticIntersector intersect = makeIntersector();
    while (out.count < ACOUSTIC_MAX_LAYERS) {
        const AcousticHitRecord hit = record(intersect.intersect(query, scene), instances, triangleFacets);
        if (hit.instance == ACOUSTIC_MISS) break;
        out.hits[out.count++] = hit;
        query.min_distance = hit.t + params.selfHitOffset;
    }
    hits[index] = out;
}

kernel void acoustic_discovery(constant AcousticTraceParams& params [[buffer(0)]],
                               const device AcousticDiscoverySeed* seeds [[buffer(1)]],
                               device AcousticDiscoveryVertex* vertices [[buffer(2)]],
                               const device AcousticInstanceRecord* instances [[buffer(3)]],
                               const device uint* triangleFacets [[buffer(4)]],
                               instance_acceleration_structure scene [[buffer(5)]],
                               const device AcousticFacetRecord* facets [[buffer(6)]],
                               uint index [[thread_position_in_grid]]) {
    const AcousticDiscoverySeed seed = seeds[index];
    uint state = seed.seed;
    float3 origin(seed.originX, seed.originY, seed.originZ);
    float3 direction(seed.directionX, seed.directionY, seed.directionZ);
    AcousticIntersector intersect = makeIntersector();
    const uint base = index * params.maxOrder;
    uint written = 0;
    while (written < params.maxOrder) {
        ray query;
        query.origin = origin;
        query.direction = direction;
        query.min_distance = params.selfHitOffset;
        query.max_distance = params.maxDistance;
        const AcousticHitRecord hit = record(intersect.intersect(query, scene), instances, triangleFacets);
        if (hit.instance == ACOUSTIC_MISS) break;
        AcousticDiscoveryVertex reached;
        reached.instance = hit.instance;
        reached.facet = hit.facet;
        reached.pad0 = reached.pad1 = 0;
        vertices[base + written++] = reached;
        const AcousticInstanceRecord instance = instances[hit.instance];
        const AcousticFacetRecord facet = facets[instance.firstFacet + hit.facet];
        const float3 object(facet.normalX, facet.normalY, facet.normalZ);
        const float3 normal = normalize(float3(dot(float3(instance.normalRows[0][0], instance.normalRows[0][1], instance.normalRows[0][2]), object),
                                               dot(float3(instance.normalRows[1][0], instance.normalRows[1][1], instance.normalRows[1][2]), object),
                                               dot(float3(instance.normalRows[2][0], instance.normalRows[2][1], instance.normalRows[2][2]), object)));
        origin = origin + direction * hit.t;
        direction = normalize(acousticNextDirection(direction, normal, facet.scattering, &state));
    }
    if (written < params.maxOrder) {
        AcousticDiscoveryVertex end;
        end.instance = ACOUSTIC_MISS;
        end.facet = 0;
        end.pad0 = end.pad1 = 0;
        vertices[base + written] = end;
    }
}
