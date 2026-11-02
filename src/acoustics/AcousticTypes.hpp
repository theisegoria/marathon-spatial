// Plain-data vocabulary of the acoustics module: octave bands, materials, transforms and the descriptions a host passes in.
#pragma once

#include <array>
#include <cmath>
#include <cstdint>
#include <simd/simd.h>
#include <span>
#include <vector>

namespace mspatial::acoustics {

inline constexpr std::size_t kBandCount = 8;
inline constexpr std::array<float, kBandCount> kBandCenterHz{62.5f, 125.f, 250.f, 500.f, 1000.f, 2000.f, 4000.f, 8000.f};
using Bands = std::array<float, kBandCount>;

inline constexpr float kDefaultSpeedOfSound = 343.f;
inline constexpr std::uint32_t kNone = 0xFFFFFFFFu;

using Vec3 = simd_float3;
inline Vec3 vec3(float x, float y, float z) { return simd_make_float3(x, y, z); }
inline float dot(Vec3 a, Vec3 b) { return simd_dot(a, b); }
inline Vec3 cross(Vec3 a, Vec3 b) { return simd_cross(a, b); }
inline float length(Vec3 a) { return simd_length(a); }
inline Vec3 normalize(Vec3 a) { return simd_normalize(a); }

inline Bands uniformBands(float value) {
    Bands out;
    out.fill(value);
    return out;
}

// Energy fractions per octave band, measured per surface crossing or per reflection.
//  absorption:   energy that does not come back. By the usual definition it includes what is transmitted.
//  scattering:   fraction of the reflected energy leaving non-specularly (ISO 17497).
//  transmission: fraction of the incident energy that gets through this one surface.
struct Material {
    Bands absorption = uniformBands(0.1f);
    Bands scattering = uniformBands(0.05f);
    Bands transmission = uniformBands(0.f);
};

// Affine object-to-world transform in row-major 3x4 form: world.x = rows[0][0..2] dot p + rows[0][3].
struct Transform {
    float rows[3][4] = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}};

    Vec3 point(Vec3 p) const {
        return vec3(rows[0][0] * p.x + rows[0][1] * p.y + rows[0][2] * p.z + rows[0][3],
                    rows[1][0] * p.x + rows[1][1] * p.y + rows[1][2] * p.z + rows[1][3],
                    rows[2][0] * p.x + rows[2][1] * p.y + rows[2][2] * p.z + rows[2][3]);
    }
    Vec3 direction(Vec3 d) const {
        return vec3(rows[0][0] * d.x + rows[0][1] * d.y + rows[0][2] * d.z,
                    rows[1][0] * d.x + rows[1][1] * d.y + rows[1][2] * d.z,
                    rows[2][0] * d.x + rows[2][1] * d.y + rows[2][2] * d.z);
    }
    static Transform translation(Vec3 t) {
        Transform out;
        out.rows[0][3] = t.x;
        out.rows[1][3] = t.y;
        out.rows[2][3] = t.z;
        return out;
    }
};

// Triangles as flat arrays. positions holds x,y,z triples; triangleMaterials holds one material id per triangle.
struct MeshDesc {
    std::span<const float> positions;
    std::span<const std::uint32_t> indices;
    std::span<const std::uint32_t> triangleMaterials;
};

struct Box {
    Vec3 min = vec3(INFINITY, INFINITY, INFINITY);
    Vec3 max = vec3(-INFINITY, -INFINITY, -INFINITY);

    void grow(Vec3 p) {
        min = simd_min(min, p);
        max = simd_max(max, p);
    }
    void grow(const Box& b) {
        min = simd_min(min, b.min);
        max = simd_max(max, b.max);
    }
    bool contains(Vec3 p, float margin = 0.f) const {
        return p.x >= min.x - margin && p.y >= min.y - margin && p.z >= min.z - margin && p.x <= max.x + margin &&
               p.y <= max.y + margin && p.z <= max.z + margin;
    }
    float area() const {
        const Vec3 e = simd_max(max - min, 0.f);
        return 2.f * (e.x * e.y + e.y * e.z + e.z * e.x);
    }
    float volume() const {
        const Vec3 e = simd_max(max - min, 0.f);
        return e.x * e.y * e.z;
    }
};

// A room is a region of late reverberation. A point belongs to the smallest room whose bounds hold it.
// With measured values at zero the simulator probes the live geometry from `probe` with rays, which follows moving doors and
// destruction; set them when the author knows better, or when the room is too far from star-shaped around one probe.
struct RoomDesc {
    std::vector<Box> bounds;
    Vec3 probe = vec3(0, 0, 0);
    float volume = 0.f;                 // m^3, 0 to measure.
    float surfaceArea = 0.f;            // m^2, 0 to measure.
    Bands absorptionArea = uniformBands(0.f);  // Equivalent absorption area in m^2 per band, used when surfaceArea is set.
};

// A rectangular opening between two rooms. normal points from room a to room b; halfU and halfV span the opening.
struct PortalDesc {
    std::uint32_t roomA = kNone;
    std::uint32_t roomB = kNone;
    Vec3 center = vec3(0, 0, 0);
    Vec3 normal = vec3(1, 0, 0);
    Vec3 halfU = vec3(0, 1, 0);
    Vec3 halfV = vec3(0, 0, 1);
    float halfThickness = 0.1f;  // Openness rays run this far to each side, so a door anywhere in the frame blocks them.

    float area() const { return 4.f * length(cross(halfU, halfV)); }
};

struct SourceDesc {
    Vec3 position = vec3(0, 0, 0);
    float radius = 0.1f;  // Physical extent of the emitter: partial occlusion is measured over this disc.
};

// ISO 9613-1 attenuation of sound in air, in energy per metre (the coefficient m in E = E0 exp(-m d)).
Bands airAttenuation(float celsius = 20.f, float relativeHumidityPercent = 50.f, float pressureKilopascal = 101.325f);

}  // namespace mspatial::acoustics
