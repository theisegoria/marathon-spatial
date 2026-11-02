// Where acoustic rays run. The simulator fills a batch of plain ray records, submits it, and picks the result up on a later tick; it
// never waits. CpuTraceBackend answers on the calling thread for tests and for machines with no GPU; MetalTraceBackend runs the
// same queries on its own Metal 4 queue.
#pragma once

#include "acoustics/AcousticScene.hpp"
#include "../../shaders/AcousticTraceShared.h"

#include <memory>

namespace mspatial::acoustics {

struct TraceCapacity {
    std::uint32_t closestRays = 2048;
    std::uint32_t layeredRays = 512;
    std::uint32_t discoveryRays = 256;
    std::uint32_t maxOrder = 4;
};

// Pointers stay valid from acquire to release. Metal batches point into shared GPU memory.
struct TraceBatch {
    std::uint32_t slot = 0;  // Which in-flight slot this is, for the caller's own per-batch records.
    bool failed = false;     // The GPU rejected or lost this batch; its hit records mean nothing.
    AcousticRayRecord* closestRays = nullptr;
    AcousticHitRecord* closestHits = nullptr;
    std::uint32_t closestCount = 0;
    AcousticRayRecord* layeredRays = nullptr;
    AcousticLayeredHit* layeredHits = nullptr;
    std::uint32_t layeredCount = 0;
    AcousticDiscoverySeed* discoverySeeds = nullptr;
    AcousticDiscoveryVertex* discoveryVertices = nullptr;  // discoveryCount * maxOrder entries.
    std::uint32_t discoveryCount = 0;
};

inline constexpr std::uint32_t kMaxBatchesInFlight = 3;

class TraceBackend {
public:
    virtual ~TraceBackend() = default;
    virtual const TraceCapacity& capacity() const = 0;
    // Takes the current instance transforms and any new meshes. Called before acquire on the simulator's thread.
    virtual void syncScene(const Scene& scene) = 0;
    // A free batch, or nullptr when every batch is in flight.
    virtual TraceBatch* acquireBatch() = 0;
    virtual void submit(TraceBatch* batch) = 0;
    // A finished batch in submission order, or nullptr. Never blocks.
    virtual TraceBatch* pollCompleted() = 0;
    virtual void release(TraceBatch* batch) = 0;
    // Rays start this far along so they do not hit the surface they left.
    virtual float selfHitOffset() const { return 1e-3f; }
};

std::unique_ptr<TraceBackend> makeCpuTraceBackend(const TraceCapacity& capacity = {});

}  // namespace mspatial::acoustics
