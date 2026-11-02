#include "acoustics/AcousticTrace.hpp"

#include <array>
#include <deque>

namespace mspatial::acoustics {
namespace {

Vec3 origin(const AcousticRayRecord& r) { return vec3(r.originX, r.originY, r.originZ); }
Vec3 direction(const AcousticRayRecord& r) { return vec3(r.directionX, r.directionY, r.directionZ); }

AcousticHitRecord toRecord(const RayHit& hit) { return {hit.instance == kNone ? ACOUSTIC_MISS : hit.instance, hit.facet, hit.t, 0.f}; }

class CpuTraceBackend final : public TraceBackend {
public:
    explicit CpuTraceBackend(const TraceCapacity& capacity) : capacity_(capacity) {
        for (std::uint32_t i = 0; i < kMaxBatchesInFlight; ++i) {
            Storage& s = storage_[i];
            s.closestRays.resize(capacity.closestRays);
            s.closestHits.resize(capacity.closestRays);
            s.layeredRays.resize(capacity.layeredRays);
            s.layeredHits.resize(capacity.layeredRays);
            s.seeds.resize(capacity.discoveryRays);
            s.vertices.resize(capacity.discoveryRays * capacity.maxOrder);
            s.batch.slot = i;
            s.batch.closestRays = s.closestRays.data();
            s.batch.closestHits = s.closestHits.data();
            s.batch.layeredRays = s.layeredRays.data();
            s.batch.layeredHits = s.layeredHits.data();
            s.batch.discoverySeeds = s.seeds.data();
            s.batch.discoveryVertices = s.vertices.data();
            free_.push_back(&s.batch);
        }
    }

    const TraceCapacity& capacity() const override { return capacity_; }

    void syncScene(const Scene& scene) override {
        scene_ = &scene;
        const std::uint64_t revision = scene.instanceRevision() + scene.meshRevision();
        if (revision == revision_) return;
        revision_ = revision;
        intersector_.rebuild(scene);
    }

    TraceBatch* acquireBatch() override {
        if (free_.empty()) return nullptr;
        TraceBatch* batch = free_.front();
        free_.pop_front();
        batch->closestCount = batch->layeredCount = batch->discoveryCount = 0;
        batch->failed = false;
        return batch;
    }

    void submit(TraceBatch* batch) override {
        for (std::uint32_t i = 0; i < batch->closestCount; ++i) {
            const AcousticRayRecord& r = batch->closestRays[i];
            batch->closestHits[i] = toRecord(intersector_.closest(*scene_, origin(r), direction(r), r.tMin, r.tMax));
        }
        for (std::uint32_t i = 0; i < batch->layeredCount; ++i) layered(batch->layeredRays[i], batch->layeredHits[i]);
        for (std::uint32_t i = 0; i < batch->discoveryCount; ++i) discovery(batch->discoverySeeds[i], batch->discoveryVertices + i * capacity_.maxOrder);
        done_.push_back(batch);
    }

    TraceBatch* pollCompleted() override {
        if (done_.empty()) return nullptr;
        TraceBatch* batch = done_.front();
        done_.pop_front();
        return batch;
    }

    void release(TraceBatch* batch) override { free_.push_back(batch); }

private:
    struct Storage {
        std::vector<AcousticRayRecord> closestRays, layeredRays;
        std::vector<AcousticHitRecord> closestHits;
        std::vector<AcousticLayeredHit> layeredHits;
        std::vector<AcousticDiscoverySeed> seeds;
        std::vector<AcousticDiscoveryVertex> vertices;
        TraceBatch batch;
    };

    void layered(const AcousticRayRecord& ray, AcousticLayeredHit& out) const {
        out = {};
        float tMin = ray.tMin;
        while (out.count < ACOUSTIC_MAX_LAYERS) {
            const RayHit hit = intersector_.closest(*scene_, origin(ray), direction(ray), tMin, ray.tMax);
            if (hit.instance == kNone) break;
            out.hits[out.count++] = toRecord(hit);
            tMin = hit.t + selfHitOffset();
        }
    }

    void discovery(const AcousticDiscoverySeed& seed, AcousticDiscoveryVertex* vertices) const {
        std::uint32_t state = seed.seed;
        Vec3 from = vec3(seed.originX, seed.originY, seed.originZ);
        Vec3 along = vec3(seed.directionX, seed.directionY, seed.directionZ);
        std::uint32_t written = 0;
        while (written < capacity_.maxOrder) {
            const RayHit hit = intersector_.closest(*scene_, from, along, selfHitOffset(), 1000.f);
            if (hit.instance == kNone) break;
            vertices[written++] = {hit.instance, hit.facet, 0.f, 0.f};
            const WorldPlane plane = scene_->worldPlane(hit.instance, hit.facet);
            const Bands& scattering = scene_->facetMaterial(hit.instance, hit.facet).scattering;
            float broadband = 0.f;
            for (float s : scattering) broadband += s / float(kBandCount);
            from = from + along * hit.t;
            along = normalize(acousticNextDirection(along, plane.normal, broadband, &state));
        }
        if (written < capacity_.maxOrder) vertices[written] = {ACOUSTIC_MISS, 0, 0.f, 0.f};
    }

    TraceCapacity capacity_;
    std::array<Storage, kMaxBatchesInFlight> storage_;
    std::deque<TraceBatch*> free_, done_;
    SceneIntersector intersector_;
    const Scene* scene_ = nullptr;
    std::uint64_t revision_ = 0;
};

}  // namespace

std::unique_ptr<TraceBackend> makeCpuTraceBackend(const TraceCapacity& capacity) { return std::make_unique<CpuTraceBackend>(capacity); }

}  // namespace mspatial::acoustics
