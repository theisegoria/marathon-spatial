// Edge diffraction around live geometry: when something stands between a source and the listener, sound still bends around its edges.
//
// Edges come from the meshes themselves: where two faces meet at an angle, or a thin surface ends. Which side of an edge is open air is
// decided per query, not from winding: the two faces split the space around the edge into two sectors, and a path can bend around the edge
// only when the source and the listener sit in the same sector and that sector is wider than half a turn. A point in the plane of one
// of the faces, where a path round the next corner of that face arrives, counts as being on the open side. That keeps a pillar's corner and
// a free-standing wall's top, and drops the inside corners of a room.
//
// Paths go over one edge or two (around both corners of a pillar), each point placed where the path is shortest. The loss per octave band
// follows Kurze and Anderson's fit to the exact wedge solution: 5 dB at the shadow boundary, rising with the Fresnel number of the detour.
// Rays decide whether a path is clear; this file only does geometry.
#pragma once

#include "acoustics/AcousticScene.hpp"

#include <array>
#include <span>

namespace mspatial::acoustics {

inline constexpr std::uint32_t kMaxDiffractionOrder = 2;
// Bends a path may hold: edge paths use at most kMaxDiffractionOrder, routes through baked probes up to this many.
inline constexpr std::uint32_t kMaxBendPoints = 4;

// Fills mesh.edges and mesh.edgeBvh from the mesh's triangles. Vertices closer than a millimetre are one vertex, so a mesh whose faces
// carry their own copies of shared corners still has its edges.
void extractDiffractionEdges(Mesh& mesh);

struct EdgeRef {
    std::uint32_t instance = kNone;
    std::uint32_t edge = 0;
    bool operator==(const EdgeRef&) const = default;
};

struct DiffractionPath {
    std::uint32_t order = 0;
    std::array<EdgeRef, kMaxBendPoints> edges{};
    std::array<Vec3, kMaxBendPoints> points{};  // Where the path bends, on the edges.
    std::array<Vec3, kMaxBendPoints> clear{};   // The same points moved 2 cm off the wedge, where validation rays turn.
    float length = 0.f;                               // Along the path, source to listener.
    Bands amplitude{};                                // Relative to free field at 1 m: wedge loss, air and the 1/r law.
    Vec3 arrival = vec3(0, 0, 0);                     // Unit vector from the listener toward the last bend.
    Vec3 emission = vec3(0, 0, 0);                    // Unit vector from the source toward the first bend.
    std::uint64_t hash() const;
};

struct DiffractionConfig {
    float maxDetour = 6.f;             // Metres of extra path. Past it even 62.5 Hz is down 17 dB more than at the boundary.
    float maxRouteDetour = 30.f;       // The same for routes through baked probes, which go further round.
    std::uint32_t candidates = 8;      // Single edges kept, shortest detour first, to pair into two-edge paths.
    float speedOfSound = kDefaultSpeedOfSound;
    Bands air = airAttenuation();
};

// Loss in dB of one wedge for a path longer than the straight line by `detour` metres (Kurze and Anderson, 1971).
float wedgeLossDb(float detour, float frequency, float speedOfSound);

// Up to out.size() paths around edges between source and listener, the loudest at 1 kHz first. Returns how many. Nothing here traces a
// ray: a path may still be blocked by something else, which the caller's rays find out.
std::uint32_t findDiffractionPaths(const Scene& scene, Vec3 source, Vec3 listener, const DiffractionConfig& config, std::span<DiffractionPath> out);

}  // namespace mspatial::acoustics
