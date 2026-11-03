// Geometric acoustics on a small ray budget. Runs on a host thread, ticks at any rate, and publishes one AcousticSnapshot per tick.
//
// Early reflections use a two-level path reservoir:
//  - Discovery walks from the listener collect facet sequences into a cache shared by every source.
//  - Each source keeps up to kMaxPathsPerSource sequences that passed image-source validation. Geometry, delay and gain of a kept path are
//    solved analytically every tick from the live transforms; rays only confirm that it is still clear. So a few rays give exact, noise-free
//    paths: temporal reuse holds the path and the rays only decide whether it lives. Scattered energy is not traced into early paths and
//    reaches the listener through the late field.
// Occlusion and transmission use a few layered rays over the source disc. Rooms and portals are measured by rays at their live positions, so a
// moving door or a blown-out wall changes all of it within a few ticks.
#pragma once

#include "acoustics/AcousticSnapshot.hpp"
#include "acoustics/AcousticTrace.hpp"
#include "acoustics/Bake.hpp"
#include "acoustics/Diffraction.hpp"
#include "acoustics/EnergyField.hpp"
#include "acoustics/RoomAcoustics.hpp"

#include <array>
#include <memory>
#include <unordered_map>

namespace mspatial::acoustics {

inline constexpr std::uint32_t kMaxReflectionOrder = 6;

struct FacetRef {
    std::uint32_t instance = kNone;
    std::uint32_t facet = 0;
    bool operator==(const FacetRef&) const = default;
};

// A reflection sequence. facets[0] is the surface nearest the source.
struct PathKey {
    std::uint32_t order = 0;
    std::array<FacetRef, kMaxReflectionOrder> facets{};
    std::uint64_t hash() const;
};

struct ImagePath {
    float length = 0.f;
    std::array<Vec3, kMaxReflectionOrder> points{};
    Vec3 arrival = vec3(0, 0, 0);
    Vec3 emission = vec3(0, 0, 0);
    Bands amplitude{};  // Relative to free field at 1 m, including surface losses, air and the 1/r law.
};

// Solves the image-source construction for a sequence against the live planes. False when the source and listener cannot reach each other
// through these planes; visibility is the caller's business.
bool solveImagePath(const Scene& scene, const PathKey& key, Vec3 source, Vec3 listener, const Bands& air, ImagePath& out);

struct SimulatorConfig {
    float speedOfSound = kDefaultSpeedOfSound;
    Bands air = airAttenuation();
    DecayFormula formula = DecayFormula::Eyring;
    std::uint32_t discoveryRays = 96;
    std::uint32_t sequenceChecksPerSource = 128;    // Analytic filter only, no rays.
    std::uint32_t newValidationsPerSource = 16;
    std::uint32_t directRays = 5;
    std::uint32_t portalRaysPerTick = 4;
    std::uint32_t probeRaysPerTick = 16;
    std::uint32_t probeWindow = 512;
    float probeMaxDistance = 40.f;
    float smoothingSeconds = 0.04f;
    float minimumPathGain = 1e-3f;
    float replacementMargin = 1.25f;  // A new path must be this much stronger in energy to retire the weakest one.
    std::uint32_t missTolerance = 2;  // Consecutive failed validations before a path fades.
    std::uint32_t sequenceCacheSize = 4096;
    // The listener's traced energy field, for convolution reverb. paths = 0 leaves it off; its rays take at most half of what the
    // direct and room rays leave, so early reflections always keep a share.
    EnergyTracerConfig energy{.paths = 0};
    // Edge diffraction for sources whose direct line is blocked. Each such source checks this many paths a tick, one ray per leg.
    bool diffraction = true;
    std::uint32_t diffractionChecksPerSource = 3;
    DiffractionConfig bends{};
    // Routes through a baked probe graph, when a bake is set: new candidates tried per blocked source each tick.
    std::uint32_t pathingChecksPerSource = 2;
};

class Simulator {
public:
    // backend must live as long as this object. Capacity checks happen here: rays beyond the backend's capacity are simply not issued.
    Simulator(Scene& scene, TraceBackend& backend, const SimulatorConfig& config = {});

    std::uint32_t addSource(const SourceDesc& desc);
    void setSource(std::uint32_t id, const SourceDesc& desc);
    void removeSource(std::uint32_t id);
    void setListener(Vec3 position) { listener_ = position; }
    Vec3 listener() const { return listener_; }
    // Routes round obstacles through this bake's probes, alongside edge diffraction. nullptr stops using one.
    void setBake(std::shared_ptr<const AcousticBake> bake) { bake_ = std::move(bake); }
    const AcousticBake* bake() const { return bake_.get(); }

    // One tick: takes finished ray results, refreshes everything analytic, publishes a snapshot, and queues the next rays.
    void tick(float dt);

    // The audio thread is the one reader of this exchange.
    SnapshotExchange<AcousticSnapshot>& snapshots() { return exchange_; }
    // The snapshot the last tick published, for the thread that calls tick. Stays intact until tick runs twice more; null before the first tick.
    const AcousticSnapshot* lastPublished() const { return published_; }
    std::uint32_t issuedClosestRays() const { return issuedClosest_; }
    // The listener's energy field, read on the thread that calls tick.
    const EnergyTracer& energyTracer() const { return energy_; }
    std::size_t knownSequences() const { return sequenceCount_; }

private:
    struct PathSlot {
        PathKey key;
        std::uint64_t id = 0;
        ImagePath path;
        Bands smoothed{};
        std::uint32_t misses = 0;
        bool visible = false;
        bool retiring = false;  // Lost its place to a stronger path: it only fades from here.
    };
    struct BendSlot {
        std::uint64_t id = 0;
        DiffractionPath path;
        Bands smoothed{};
        std::uint32_t misses = 0;
        bool visible = false;
    };
    struct SourceState {
        SourceDesc desc;
        bool active = false;
        std::uint32_t generation = 0;
        std::uint32_t room = kNone;
        std::uint32_t cursor = 0;
        std::uint32_t confirmCursor = 0;
        Bands directTarget = uniformBands(1.f);
        Bands directEnergy = uniformBands(1.f);
        float visibilityTarget = 1.f;
        float visibility = 1.f;
        std::array<PathSlot, kMaxPathsPerSource> slots;
        std::array<BendSlot, kMaxDiffractionPaths> bends;
    };
    struct BendJob {
        std::uint32_t source, generation, firstRay, rays;
        DiffractionPath path;
        std::uint32_t slot;  // The slot it confirms, or kNone for a new candidate.
    };
    struct PathJob {
        std::uint32_t source, generation, firstRay;
        PathKey key;
        std::array<float, kMaxReflectionOrder + 1> lengths;
        bool candidate;
    };
    struct DirectJob {
        std::uint32_t source, generation, firstRay, rays;
    };
    struct SampleJob {
        std::uint32_t owner, sample, ray;
    };
    struct Jobs {
        std::vector<PathJob> paths;
        std::vector<DirectJob> direct;
        std::vector<SampleJob> portals, probes;
        std::vector<Vec3> probeDirections;
        std::vector<EnergyJob> energy;
        std::vector<BendJob> bends;
        std::uint32_t discovery = 0;
    };
    struct ProbeSample {
        float cubedRange = 0.f;
        float area = 0.f;
        Bands absorption{};
        bool filled = false;
    };
    struct RoomState {
        std::vector<ProbeSample> samples;
        std::uint32_t cursor = 0;
    };
    struct PortalState {
        std::array<Bands, 16> samples;
        std::array<bool, 16> filled{};
        std::uint32_t cursor = 0;
        Bands open = uniformBands(1.f);
    };

    void prepare(TraceBatch& batch);
    void consume(const TraceBatch& batch);
    void consumePaths(const TraceBatch& batch, const Jobs& jobs);
    void consumeDirect(const TraceBatch& batch, const Jobs& jobs);
    void prepareBends(TraceBatch& batch, Jobs& jobs);
    void consumeBends(const TraceBatch& batch, const Jobs& jobs);
    void consumeSamples(const TraceBatch& batch, const Jobs& jobs);
    void consumeDiscovery(const TraceBatch& batch, const Jobs& jobs);
    void rememberSequence(const PathKey& key);
    bool pushPathJob(TraceBatch& batch, Jobs& jobs, std::uint32_t source, const PathKey& key, const ImagePath& path, bool candidate);
    void admit(SourceState& source, const PathKey& key, const ImagePath& path);
    void updateRooms(AcousticSnapshot& out);
    void updatePortalPaths(std::uint32_t source, SourceAcoustics& out);
    void publish(float dt);

    Scene& scene_;
    TraceBackend& backend_;
    SimulatorConfig config_;
    Vec3 listener_ = vec3(0, 0, 0);
    std::uint64_t tick_ = 0;
    std::uint32_t issuedClosest_ = 0;
    std::uint32_t fairSource_ = 0;
    std::array<SourceState, kMaxSources> sources_;
    std::array<Jobs, kMaxBatchesInFlight> jobs_;
    std::vector<RoomState> roomStates_;
    std::vector<PortalState> portalStates_;
    std::vector<PathKey> sequences_;
    std::unordered_map<std::uint64_t, std::uint32_t> sequenceIndex_;
    std::uint32_t sequenceHead_ = 0;
    std::size_t sequenceCount_ = 0;
    std::vector<float> portalDistance_;
    std::vector<std::uint32_t> portalPrevious_;
    std::vector<std::uint8_t> portalDone_;
    EnergyTracer energy_;
    std::shared_ptr<const AcousticBake> bake_;
    SnapshotExchange<AcousticSnapshot> exchange_;
    const AcousticSnapshot* published_ = nullptr;
};

}  // namespace mspatial::acoustics
