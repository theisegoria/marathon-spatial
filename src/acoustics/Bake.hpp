// Baked propagation for the static part of a scene: probes placed through the rooms, each holding the energy field a listener there would
// measure, and a graph of which probes see each other with the shortest route between every pair.
//
// The live simulator needs neither: it traces everything against the moving scene. A bake buys two things it cannot get cheaply.
//  - Reverb at once. A listener entering a room starts with the converged field of the nearest probes instead of a few hundred paths, and
//    the live field refines it as evidence comes in: the bake is prior evidence worth a fixed number of paths, so the live field takes over
//    wherever the scene has changed since.
//  - Routes of any length. Edge diffraction finds paths over one or two edges. Sound reaching a listener round three corners of a corridor
//    takes a route through the probe graph, pulled tight, and the live rays still decide whether each leg is clear.
//
// Probes come from a grid over each room, kept where they stand clear of every surface, or from the host. Grid probes that see no other probe
// (inside a solid pillar, say) are dropped. A bake is plain data and saves to a versioned little-endian file.
#pragma once

#include "acoustics/Diffraction.hpp"
#include "acoustics/EnergyField.hpp"

#include <functional>
#include <string>

namespace mspatial::acoustics {

// EdgeRef::instance of a bend at a baked probe; EdgeRef::edge is the probe.
inline constexpr std::uint32_t kBakeProbeTag = 0xFFFFFFFEu;
inline constexpr std::uint32_t kMaxBakeProbes = 4096;

struct BakeConfig {
    float spacing = 2.f;               // Metres between grid probes.
    float clearance = 0.3f;            // A grid probe this close to any surface is dropped.
    std::vector<Vec3> probes;          // Host-placed probes. When set, no grid is placed and none are dropped.
    bool reverb = true;
    bool pathing = true;
    std::uint32_t pathsPerProbe = 4096;  // Energy paths traced to the end of the field from each probe.
    // Length, bins, receiver and air of the field. paths is how many run at once on each thread; age is never applied.
    EnergyTracerConfig energy{.paths = 512};
    std::uint32_t threads = 4;
};

struct BakeView {
    std::uint32_t probes = 0;
    std::uint32_t droppedProbes = 0;   // Grid probes too near a surface or seeing no other probe.
    std::uint64_t visiblePairs = 0;
    std::uint64_t segments = 0;        // Energy path segments traced.
};

class AcousticBake {
public:
    std::vector<Vec3> positions;
    std::vector<std::uint32_t> rooms;  // Room of each probe, kNone outside every room.
    std::uint64_t sceneHash = 0;       // staticSceneHash of the scene it was baked from.
    BakeView stats;

    // Reverb. A probe's field is stored as its total per band and bin, and the share each direction takes in two windows: before and after
    // earlySeconds. The renderer already spreads direction slowly over time, so this keeps what it uses at a tenth of the size.
    std::uint32_t bins = 0;
    float binSeconds = 0.01f;
    float earlySeconds = 0.08f;
    double pathsPerProbe = 0.0;
    std::vector<float> totals;  // [probe][band][bin]
    std::vector<float> shares;  // [probe][window][direction][band]

    // Pathing.
    std::vector<std::uint64_t> visibility;  // Row a holds bit b when probes a and b see each other.
    std::vector<float> distances;           // [a][b] shortest route through probes, INFINITY when there is none.
    std::vector<std::uint16_t> next;        // [a][b] the probe after a on that route.

    std::uint32_t probeCount() const { return std::uint32_t(positions.size()); }
    bool hasReverb() const { return bins > 0 && !totals.empty(); }
    bool hasPathing() const { return !distances.empty(); }
    bool sees(std::uint32_t a, std::uint32_t b) const;

    // The field a listener at point would measure, blended from up to four of the nearest probes in that room, laid out in the given bins.
    // False when the room has no probe.
    bool fieldAt(Vec3 point, std::uint32_t room, std::uint32_t bins, float binSeconds, EnergyField& out) const;

    // Routes from source to listener through probes, at most kMaxBendPoints bends each, pulled tight over the baked visibility. Nothing here
    // traces a ray: the first and last legs are guesses until the caller's rays confirm them. Successive values of rotate walk through the
    // candidates, so a guess that fails gives way to the next. Returns how many it wrote.
    std::uint32_t findPaths(Vec3 source, Vec3 listener, const DiffractionConfig& config, std::uint32_t rotate, std::span<DiffractionPath> out) const;

    bool save(const std::string& path, std::string* error = nullptr) const;
    static bool load(const std::string& path, AcousticBake& out, std::string* error = nullptr);
};

// True when the path bends at baked probes rather than at edges.
inline bool isProbePath(const DiffractionPath& path) { return path.order > 0 && path.edges[0].instance == kBakeProbeTag; }

// Recomputes length, directions and per-band amplitude of a route through probes for the live source and listener. The route loses what one
// wedge with the route's whole detour would (Kurze and Anderson), with the 5 dB at the shadow boundary ramped in as its sharpest bend turns
// through the first 20 degrees. Probes stand in open air away from the corners, so per-bend losses would overstate a route that is only an
// estimate; this is the gentler, smooth estimate, and the edge paths keep the exact per-wedge loss.
void finishProbePath(DiffractionPath& path, Vec3 source, Vec3 listener, const DiffractionConfig& config);

// A hash of what rays see in the scene: enabled instances, their transforms, triangles and materials. Equal hashes mean a bake still fits.
std::uint64_t staticSceneHash(const Scene& scene);

// Bakes the scene. Slow: seconds to minutes, on config.threads threads, the calling one among them. progress, when set, is called on the
// calling thread with the share done, from 0 to 1.
bool bakeScene(const Scene& scene, const BakeConfig& config, AcousticBake& out, std::string* error = nullptr,
               const std::function<void(float)>& progress = {});

// Mixes a baked field into a live one as prior evidence worth priorPaths paths. Both use the live field's bins.
void blendFields(const EnergyField& live, const EnergyField& baked, double priorPaths, EnergyField& out);

}  // namespace mspatial::acoustics
