// Layouts and the discovery bounce rule shared by the Metal kernels and the CPU backend, so the two trace the same paths.
// Every record uses scalar members: no float3 padding differs between MSL and C++.
#pragma once

#ifdef __METAL_VERSION__
#include <metal_stdlib>
using namespace metal;
typedef uint acoustic_uint;
typedef float3 acoustic_float3;
#define ACOUSTIC_INLINE inline
#define ACOUSTIC_THREAD thread
#else
#include <cstdint>
#include <simd/simd.h>
typedef std::uint32_t acoustic_uint;
typedef simd_float3 acoustic_float3;
#define ACOUSTIC_INLINE inline
#define ACOUSTIC_THREAD
#endif

#define ACOUSTIC_MAX_LAYERS 4u
#define ACOUSTIC_MISS 0xFFFFFFFFu

struct AcousticRayRecord {
    float originX, originY, originZ, tMin;
    float directionX, directionY, directionZ, tMax;
};

struct AcousticHitRecord {
    acoustic_uint instance;  // Scene instance slot, ACOUSTIC_MISS for no hit.
    acoustic_uint facet;     // Facet index inside the instance's mesh.
    float t;
    float pad;
};

struct AcousticLayeredHit {
    acoustic_uint count;
    acoustic_uint pad0, pad1, pad2;
    struct AcousticHitRecord hits[ACOUSTIC_MAX_LAYERS];
};

struct AcousticDiscoverySeed {
    float originX, originY, originZ;
    acoustic_uint seed;
    float directionX, directionY, directionZ;
    acoustic_uint pad;
};

struct AcousticDiscoveryVertex {
    acoustic_uint instance;  // ACOUSTIC_MISS ends the chain.
    acoustic_uint facet;
    float pad0, pad1;
};

// One scene instance as the GPU sees it. normalRows is the inverse transpose of the linear part, which keeps reflections
// correct on scaled instances.
struct AcousticInstanceRecord {
    float normalRows[3][4];
    acoustic_uint firstTriangle;
    acoustic_uint firstFacet;
    acoustic_uint pad0, pad1;
};

// Object-space facet normal and a broadband scattering value for the discovery walk.
struct AcousticFacetRecord {
    float normalX, normalY, normalZ;
    float scattering;
};

struct AcousticTraceParams {
    acoustic_uint rayCount;
    acoustic_uint maxOrder;
    float selfHitOffset;
    float maxDistance;
};

ACOUSTIC_INLINE acoustic_uint acousticRandom(ACOUSTIC_THREAD acoustic_uint* state) {
    const acoustic_uint s = *state * 747796405u + 2891336453u;
    *state = s;
    const acoustic_uint word = ((s >> ((s >> 28u) + 4u)) ^ s) * 277803737u;
    return (word >> 22u) ^ word;
}

ACOUSTIC_INLINE float acousticRandomUnit(ACOUSTIC_THREAD acoustic_uint* state) {
    return float(acousticRandom(state) >> 8u) * (1.0f / 16777216.0f);
}

ACOUSTIC_INLINE float acousticDot(acoustic_float3 a, acoustic_float3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }

// The next direction of a discovery ray that hit a surface of normal `normal`. A hit reflects specularly, or with probability
// `scattering` leaves in a cosine lobe, which reaches facet sequences a purely specular walk would never find.
ACOUSTIC_INLINE acoustic_float3 acousticNextDirection(acoustic_float3 direction, acoustic_float3 normal, float scattering,
                                                      ACOUSTIC_THREAD acoustic_uint* state) {
    const float incidence = acousticDot(direction, normal);
    const acoustic_float3 facing = incidence > 0.0f ? -normal : normal;
    if (acousticRandomUnit(state) >= scattering) return direction - 2.0f * incidence * normal;
    const float u = acousticRandomUnit(state);
    const float v = acousticRandomUnit(state) * 6.28318530718f;
    const float r = sqrt(u);
#ifdef __METAL_VERSION__
    const float cosV = cos(v), sinV = sin(v), height = sqrt(max(0.0f, 1.0f - u));
#else
    const float cosV = cosf(v), sinV = sinf(v), height = sqrtf(fmaxf(0.0f, 1.0f - u));
#endif
    const float sign = facing.z >= 0.0f ? 1.0f : -1.0f;
    const float a = -1.0f / (sign + facing.z);
    const float b = facing.x * facing.y * a;
    acoustic_float3 tangent, bitangent;
    tangent.x = 1.0f + sign * facing.x * facing.x * a;
    tangent.y = sign * b;
    tangent.z = -sign * facing.x;
    bitangent.x = b;
    bitangent.y = sign + facing.y * facing.y * a;
    bitangent.z = -facing.y;
    return r * cosV * tangent + r * sinV * bitangent + height * facing;
}
