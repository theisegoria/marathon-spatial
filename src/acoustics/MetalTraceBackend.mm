// Runs acoustic ray batches on a private MTL4CommandQueue. The batch memory is shared storage: the simulator writes rays into it directly and
// reads hits straight back. Completion is observed by polling a shared event, so nothing here ever blocks, except teardown.
//
// Instance acceleration structures are rebuilt for every batch that follows a transform change, using the indirect instance form, with
// PreferFastBuild: the instance count is small, and it keeps moving doors exact. Bottom-level structures only change when meshes change.
#import <Metal/Metal.h>

#include "acoustics/MetalTraceBackend.hpp"

#include <atomic>
#include <deque>
#include <memory>
#include <unordered_map>

#if !__has_feature(objc_arc)
#error "MetalTraceBackend.mm must be compiled with -fobjc-arc"
#endif

namespace mspatial::acoustics {
namespace {

constexpr NSUInteger kBindParams = 0, kBindRays = 1, kBindHits = 2, kBindInstances = 3, kBindTriangleFacets = 4, kBindScene = 5, kBindFacets = 6;
constexpr MTLResourceOptions kShared = MTLResourceStorageModeShared | MTLResourceHazardTrackingModeUntracked;
constexpr float kDiscoveryRange = 1000.f;  // Same as the CPU backend.

struct MeshEntry {
    std::uint64_t serial = 0;
    id<MTLBuffer> vertices;
    id<MTLBuffer> indices;
    id<MTLAccelerationStructure> structure;
    MTL4PrimitiveAccelerationStructureDescriptor* descriptor = nil;
    NSUInteger scratchSize = 0;
    bool built = false;
};

// Anything a batch in flight may still use is held here until the submission counter passes it.
struct Retired {
    std::uint64_t safeAfter = 0;
    std::vector<id<MTLAllocation>> allocations;
};

struct BatchStorage {
    TraceBatch batch;
    id<MTLBuffer> closestRays, closestHits, layeredRays, layeredHits, seeds, vertices, params;
    id<MTLBuffer> instanceRecords, instanceDescriptors, scratch;
    id<MTLAccelerationStructure> structure;
    id<MTL4CommandAllocator> allocator;
    id<MTL4CommandBuffer> commands;
    id<MTL4ArgumentTable> closestTable, layeredTable, discoveryTable;
    std::uint64_t completion = 0;
    std::uint64_t builtRevision = 0;
    NSUInteger instanceCount = 0;
    std::shared_ptr<std::atomic<bool>> failed = std::make_shared<std::atomic<bool>>(false);
};

class MetalTraceBackend final : public TraceBackend {
public:
    MetalTraceBackend(const MetalTraceConfig& config, std::string* error) : capacity_(config.capacity) {
        device_ = config.device ? (__bridge id<MTLDevice>)config.device : MTLCreateSystemDefaultDevice();
        const auto fail = [&](NSString* message) {
            if (error) *error = message.UTF8String;
        };
        if (!device_ || !device_.supportsRaytracing || ![device_ supportsFamily:MTLGPUFamilyMetal4] || ![device_ supportsFamily:MTLGPUFamilyApple9]) {
            fail(@"Metal 4 ray tracing on Apple family 9 or later is required");
            return;
        }
        NSError* nsError = nil;
        id<MTLLibrary> library = [device_ newLibraryWithURL:[NSURL fileURLWithPath:@(config.libraryPath.c_str())] error:&nsError];
        if (!library) {
            fail(nsError.localizedDescription ?: @"could not load the acoustic trace library");
            return;
        }
        MTL4CompilerDescriptor* compilerDescriptor = [MTL4CompilerDescriptor new];
        compilerDescriptor.label = @"Acoustics compiler";
        id<MTL4Compiler> compiler = [device_ newCompilerWithDescriptor:compilerDescriptor error:&nsError];
        if (!compiler) {
            fail(nsError.localizedDescription ?: @"no Metal 4 compiler");
            return;
        }
        const auto pipeline = [&](NSString* name) -> id<MTLComputePipelineState> {
            MTL4LibraryFunctionDescriptor* function = [MTL4LibraryFunctionDescriptor new];
            function.library = library;
            function.name = name;
            MTL4ComputePipelineDescriptor* descriptor = [MTL4ComputePipelineDescriptor new];
            descriptor.computeFunctionDescriptor = function;
            descriptor.label = name;
            return [compiler newComputePipelineStateWithDescriptor:descriptor compilerTaskOptions:nil error:&nsError];
        };
        closest_ = pipeline(@"acoustic_closest");
        layered_ = pipeline(@"acoustic_layered");
        discovery_ = pipeline(@"acoustic_discovery");
        if (!closest_ || !layered_ || !discovery_) {
            fail(nsError.localizedDescription ?: @"could not build the acoustic pipelines");
            return;
        }
        queue_ = [device_ newMTL4CommandQueue];
        event_ = [device_ newSharedEvent];
        MTLResidencySetDescriptor* residency = [MTLResidencySetDescriptor new];
        residency.label = @"Acoustics";
        residency.initialCapacity = 64;
        residency_ = [device_ newResidencySetWithDescriptor:residency error:&nsError];
        if (!queue_ || !event_ || !residency_) {
            fail(nsError.localizedDescription ?: @"could not create the acoustics queue");
            return;
        }
        [queue_ addResidencySet:residency_];
        for (std::uint32_t i = 0; i < kMaxBatchesInFlight; ++i) {
            if (!makeBatch(storage_[i], i)) {
                fail(@"could not allocate acoustic batch memory");
                return;
            }
            free_.push_back(&storage_[i]);
        }
        [residency_ commit];
        [residency_ requestResidency];
        ready_ = true;
    }

    ~MetalTraceBackend() override {
        // Command buffers do not retain resources, so teardown is the one place this object waits. It runs on the host thread, never on audio.
        if (event_ && submitted_ > 0) [event_ waitUntilSignaledValue:submitted_ timeoutMS:5000];
        if (residency_) [residency_ endResidency];
    }

    bool ready() const { return ready_; }
    const TraceCapacity& capacity() const override { return capacity_; }

    void syncScene(const Scene& scene) override {
        scene_ = &scene;
        if (scene.meshRevision() != meshRevision_ || scene.materialRevision() != materialRevision_) {
            meshRevision_ = scene.meshRevision();
            materialRevision_ = scene.materialRevision();
            syncMeshes(scene);
            syncTables(scene);
            [residency_ commit];
        }
        collectRetired();
    }

    TraceBatch* acquireBatch() override {
        if (free_.empty()) return nullptr;
        BatchStorage* storage = free_.front();
        free_.pop_front();
        storage->batch.closestCount = storage->batch.layeredCount = storage->batch.discoveryCount = 0;
        storage->batch.failed = false;
        return &storage->batch;
    }

    void submit(TraceBatch* batch) override {
        BatchStorage& s = storage_[batch->slot];
        s.failed->store(false);
        if (!encode(s)) {
            // Nothing to trace against: no instances yet. All rays miss, which is what an empty world does.
            markMisses(s);
            s.completion = 0;
            done_.push_back(&s);
            return;
        }
        s.completion = ++submitted_;
        MTL4CommitOptions* options = [MTL4CommitOptions new];
        std::shared_ptr<std::atomic<bool>> failed = s.failed;
        [options addFeedbackHandler:^(id<MTL4CommitFeedback> feedback) {
            if (feedback.error) failed->store(true);
        }];
        id<MTL4CommandBuffer> buffers[] = {s.commands};
        [queue_ commit:buffers count:1 options:options];
        [queue_ signalEvent:event_ value:s.completion];
        done_.push_back(&s);
    }

    TraceBatch* pollCompleted() override {
        if (done_.empty()) return nullptr;
        BatchStorage* s = done_.front();
        if (s->completion != 0 && event_.signaledValue < s->completion) return nullptr;
        done_.pop_front();
        s->batch.failed = s->failed->load();
        return &s->batch;
    }

    void release(TraceBatch* batch) override { free_.push_back(&storage_[batch->slot]); }

private:
    id<MTLBuffer> shared(NSUInteger bytes) {
        id<MTLBuffer> buffer = [device_ newBufferWithLength:std::max<NSUInteger>(bytes, 16) options:kShared];
        if (buffer) [residency_ addAllocation:buffer];
        return buffer;
    }

    id<MTL4ArgumentTable> table() {
        MTL4ArgumentTableDescriptor* descriptor = [MTL4ArgumentTableDescriptor new];
        descriptor.maxBufferBindCount = kBindFacets + 1;
        descriptor.initializeBindings = YES;
        NSError* error = nil;
        return [device_ newArgumentTableWithDescriptor:descriptor error:&error];
    }

    bool makeBatch(BatchStorage& s, std::uint32_t slot) {
        s.closestRays = shared(capacity_.closestRays * sizeof(AcousticRayRecord));
        s.closestHits = shared(capacity_.closestRays * sizeof(AcousticHitRecord));
        s.layeredRays = shared(capacity_.layeredRays * sizeof(AcousticRayRecord));
        s.layeredHits = shared(capacity_.layeredRays * sizeof(AcousticLayeredHit));
        s.seeds = shared(capacity_.discoveryRays * sizeof(AcousticDiscoverySeed));
        s.vertices = shared(capacity_.discoveryRays * capacity_.maxOrder * sizeof(AcousticDiscoveryVertex));
        s.params = shared(sizeof(AcousticTraceParams));
        s.allocator = [device_ newCommandAllocator];
        s.commands = [device_ newCommandBuffer];
        s.closestTable = table();
        s.layeredTable = table();
        s.discoveryTable = table();
        if (!s.closestRays || !s.closestHits || !s.layeredRays || !s.layeredHits || !s.seeds || !s.vertices || !s.params || !s.allocator || !s.commands ||
            !s.closestTable || !s.layeredTable || !s.discoveryTable)
            return false;
        s.batch.slot = slot;
        s.batch.closestRays = static_cast<AcousticRayRecord*>(s.closestRays.contents);
        s.batch.closestHits = static_cast<AcousticHitRecord*>(s.closestHits.contents);
        s.batch.layeredRays = static_cast<AcousticRayRecord*>(s.layeredRays.contents);
        s.batch.layeredHits = static_cast<AcousticLayeredHit*>(s.layeredHits.contents);
        s.batch.discoverySeeds = static_cast<AcousticDiscoverySeed*>(s.seeds.contents);
        s.batch.discoveryVertices = static_cast<AcousticDiscoveryVertex*>(s.vertices.contents);
        return true;
    }

    void retire(std::vector<id<MTLAllocation>> allocations) {
        retired_.push_back({submitted_, std::move(allocations)});
    }

    void collectRetired() {
        const std::uint64_t finished = event_.signaledValue;
        while (!retired_.empty() && retired_.front().safeAfter <= finished) {
            for (id<MTLAllocation> allocation : retired_.front().allocations) [residency_ removeAllocation:allocation];
            retired_.pop_front();
        }
    }

    void syncMeshes(const Scene& scene) {
        const auto& meshes = scene.meshes();
        for (auto it = entries_.begin(); it != entries_.end();) {
            const std::uint32_t id = it->first;
            if (id < meshes.size() && meshes[id].live && meshes[id].serial == it->second.serial) {
                ++it;
                continue;
            }
            retire({it->second.vertices, it->second.indices, it->second.structure});
            it = entries_.erase(it);
        }
        for (std::uint32_t id = 0; id < meshes.size(); ++id) {
            const Mesh& mesh = meshes[id];
            if (!mesh.live || mesh.indices.empty() || entries_.contains(id)) continue;
            MeshEntry entry;
            entry.serial = mesh.serial;
            entry.vertices = shared(mesh.vertices.size() * sizeof(float) * 3);
            entry.indices = shared(mesh.indices.size() * sizeof(std::uint32_t));
            auto* packed = static_cast<float*>(entry.vertices.contents);
            for (std::size_t v = 0; v < mesh.vertices.size(); ++v) {
                packed[v * 3] = mesh.vertices[v].x;
                packed[v * 3 + 1] = mesh.vertices[v].y;
                packed[v * 3 + 2] = mesh.vertices[v].z;
            }
            std::memcpy(entry.indices.contents, mesh.indices.data(), mesh.indices.size() * sizeof(std::uint32_t));
            MTL4AccelerationStructureTriangleGeometryDescriptor* geometry = [MTL4AccelerationStructureTriangleGeometryDescriptor new];
            geometry.vertexBuffer = MTL4BufferRangeMake(entry.vertices.gpuAddress, entry.vertices.length);
            geometry.vertexFormat = MTLAttributeFormatFloat3;
            geometry.vertexStride = sizeof(float) * 3;
            geometry.indexBuffer = MTL4BufferRangeMake(entry.indices.gpuAddress, entry.indices.length);
            geometry.indexType = MTLIndexTypeUInt32;
            geometry.triangleCount = mesh.indices.size() / 3;
            geometry.opaque = YES;
            entry.descriptor = [MTL4PrimitiveAccelerationStructureDescriptor new];
            entry.descriptor.geometryDescriptors = @[geometry];
            const MTLAccelerationStructureSizes sizes = [device_ accelerationStructureSizesWithDescriptor:entry.descriptor];
            entry.structure = [device_ newAccelerationStructureWithSize:sizes.accelerationStructureSize];
            entry.scratchSize = sizes.buildScratchBufferSize;
            [residency_ addAllocation:entry.structure];
            entries_.emplace(id, std::move(entry));
        }
    }

    void syncTables(const Scene& scene) {
        std::size_t triangles = 0, facets = 0;
        firstTriangle_.assign(scene.meshes().size(), 0);
        firstFacet_.assign(scene.meshes().size(), 0);
        for (std::size_t m = 0; m < scene.meshes().size(); ++m) {
            firstTriangle_[m] = std::uint32_t(triangles);
            firstFacet_[m] = std::uint32_t(facets);
            triangles += scene.meshes()[m].triangleFacet.size();
            facets += scene.meshes()[m].facets.size();
        }
        std::vector<id<MTLAllocation>> old;
        if (triangleFacets_) old.push_back(triangleFacets_);
        if (facets_) old.push_back(facets_);
        if (!old.empty()) retire(std::move(old));
        triangleFacets_ = shared(triangles * sizeof(std::uint32_t));
        facets_ = shared(facets * sizeof(AcousticFacetRecord));
        auto* triangleOut = static_cast<std::uint32_t*>(triangleFacets_.contents);
        auto* facetOut = static_cast<AcousticFacetRecord*>(facets_.contents);
        for (std::size_t m = 0; m < scene.meshes().size(); ++m) {
            const Mesh& mesh = scene.meshes()[m];
            std::memcpy(triangleOut + firstTriangle_[m], mesh.triangleFacet.data(), mesh.triangleFacet.size() * sizeof(std::uint32_t));
            for (std::size_t f = 0; f < mesh.facets.size(); ++f) {
                float scattering = 0.f;
                for (float value : scene.materials()[mesh.facets[f].material].scattering) scattering += value / float(kBandCount);
                facetOut[firstFacet_[m] + f] = {mesh.facets[f].normal.x, mesh.facets[f].normal.y, mesh.facets[f].normal.z, scattering};
            }
        }
    }

    bool writeInstances(BatchStorage& s) {
        const auto& instances = scene_->instances();
        const NSUInteger recordBytes = instances.size() * sizeof(AcousticInstanceRecord);
        const NSUInteger descriptorBytes = instances.size() * sizeof(MTLIndirectAccelerationStructureInstanceDescriptor);
        if (s.instanceRecords.length < recordBytes) {
            if (s.instanceRecords) retire({s.instanceRecords, s.instanceDescriptors});
            s.instanceRecords = shared(recordBytes * 2);
            s.instanceDescriptors = shared(descriptorBytes * 2);
        }
        if (!s.instanceRecords || !s.instanceDescriptors) return false;
        auto* records = static_cast<AcousticInstanceRecord*>(s.instanceRecords.contents);
        auto* descriptors = static_cast<MTLIndirectAccelerationStructureInstanceDescriptor*>(s.instanceDescriptors.contents);
        s.instanceCount = 0;
        for (std::uint32_t i = 0; i < instances.size(); ++i) {
            const Instance& instance = instances[i];
            AcousticInstanceRecord record{};
            if (!instance.live || !instance.enabled) {
                records[i] = record;
                continue;
            }
            const auto entry = entries_.find(instance.mesh);
            if (entry == entries_.end()) continue;
            std::memcpy(record.normalRows, instance.normalToWorld.rows, sizeof(record.normalRows));
            record.firstTriangle = firstTriangle_[instance.mesh];
            record.firstFacet = firstFacet_[instance.mesh];
            records[i] = record;
            MTLIndirectAccelerationStructureInstanceDescriptor descriptor{};
            for (int column = 0; column < 4; ++column) {
                descriptor.transformationMatrix.columns[column].x = instance.objectToWorld.rows[0][column];
                descriptor.transformationMatrix.columns[column].y = instance.objectToWorld.rows[1][column];
                descriptor.transformationMatrix.columns[column].z = instance.objectToWorld.rows[2][column];
            }
            descriptor.options = MTLAccelerationStructureInstanceOptionOpaque | MTLAccelerationStructureInstanceOptionDisableTriangleCulling;
            descriptor.mask = 0xFFFFFFFFu;
            descriptor.userID = i;
            descriptor.accelerationStructureID = entry->second.structure.gpuResourceID;
            descriptors[s.instanceCount++] = descriptor;
        }
        return true;
    }

    void markMisses(BatchStorage& s) {
        for (std::uint32_t i = 0; i < s.batch.closestCount; ++i) s.batch.closestHits[i] = {ACOUSTIC_MISS, 0, 0.f, 0.f};
        for (std::uint32_t i = 0; i < s.batch.layeredCount; ++i) s.batch.layeredHits[i] = {};
        for (std::uint32_t i = 0; i < s.batch.discoveryCount * capacity_.maxOrder; ++i) s.batch.discoveryVertices[i] = {ACOUSTIC_MISS, 0, 0.f, 0.f};
    }

    // Records a full batch: new bottom levels, this batch's instance structure, then the three kernels.
    bool encode(BatchStorage& s) {
        if (!scene_ || !writeInstances(s) || s.instanceCount == 0) return false;

        MTL4InstanceAccelerationStructureDescriptor* top = [MTL4InstanceAccelerationStructureDescriptor new];
        top.instanceDescriptorBuffer = MTL4BufferRangeMake(s.instanceDescriptors.gpuAddress, s.instanceDescriptors.length);
        top.instanceDescriptorStride = sizeof(MTLIndirectAccelerationStructureInstanceDescriptor);
        top.instanceCount = s.instanceCount;
        top.instanceDescriptorType = MTLAccelerationStructureInstanceDescriptorTypeIndirect;
        top.usage = MTLAccelerationStructureUsagePreferFastBuild;
        const MTLAccelerationStructureSizes topSizes = [device_ accelerationStructureSizesWithDescriptor:top];

        NSUInteger scratchBytes = topSizes.buildScratchBufferSize;
        for (const auto& [id, entry] : entries_)
            if (!entry.built) scratchBytes = std::max(scratchBytes, entry.scratchSize);
        if (s.scratch.length < scratchBytes) {
            if (s.scratch) retire({s.scratch});
            s.scratch = shared(scratchBytes * 2);
        }
        if (s.structure == nil || s.structure.size < topSizes.accelerationStructureSize) {
            if (s.structure) retire({s.structure});
            s.structure = [device_ newAccelerationStructureWithSize:topSizes.accelerationStructureSize * 2];
            [residency_ addAllocation:s.structure];
        }
        if (!s.scratch || !s.structure) return false;
        [residency_ commit];

        [s.allocator reset];
        [s.commands beginCommandBufferWithAllocator:s.allocator];
        id<MTL4ComputeCommandEncoder> encoder = [s.commands computeCommandEncoder];
        const MTL4BufferRange scratch = MTL4BufferRangeMake(s.scratch.gpuAddress, s.scratch.length);
        bool builtAny = false;
        for (auto& [id, entry] : entries_) {
            if (entry.built) continue;
            if (builtAny)
                [encoder barrierAfterEncoderStages:MTLStageAccelerationStructure beforeEncoderStages:MTLStageAccelerationStructure visibilityOptions:MTL4VisibilityOptionNone];
            [encoder buildAccelerationStructure:entry.structure descriptor:entry.descriptor scratchBuffer:scratch];
            entry.built = true;
            entry.descriptor = nil;
            builtAny = true;
        }
        if (builtAny)
            [encoder barrierAfterEncoderStages:MTLStageAccelerationStructure beforeEncoderStages:MTLStageAccelerationStructure visibilityOptions:MTL4VisibilityOptionNone];
        [encoder buildAccelerationStructure:s.structure descriptor:top scratchBuffer:scratch];
        [encoder barrierAfterEncoderStages:MTLStageAccelerationStructure beforeEncoderStages:MTLStageDispatch visibilityOptions:MTL4VisibilityOptionDevice];

        auto* params = static_cast<AcousticTraceParams*>(s.params.contents);
        *params = {0, capacity_.maxOrder, selfHitOffset(), kDiscoveryRange};
        const auto bind = [&](id<MTL4ArgumentTable> table, id<MTLBuffer> rays, id<MTLBuffer> hits) {
            [table setAddress:s.params.gpuAddress atIndex:kBindParams];
            [table setAddress:rays.gpuAddress atIndex:kBindRays];
            [table setAddress:hits.gpuAddress atIndex:kBindHits];
            [table setAddress:s.instanceRecords.gpuAddress atIndex:kBindInstances];
            [table setAddress:triangleFacets_.gpuAddress atIndex:kBindTriangleFacets];
            [table setResource:s.structure.gpuResourceID atBufferIndex:kBindScene];
            [table setAddress:facets_.gpuAddress atIndex:kBindFacets];
        };
        const auto dispatch = [&](id<MTLComputePipelineState> pipeline, id<MTL4ArgumentTable> table, std::uint32_t count) {
            if (count == 0) return;
            [encoder setComputePipelineState:pipeline];
            [encoder setArgumentTable:table];
            [encoder dispatchThreads:MTLSizeMake(count, 1, 1) threadsPerThreadgroup:MTLSizeMake(std::min<NSUInteger>(pipeline.threadExecutionWidth, count), 1, 1)];
        };
        bind(s.closestTable, s.closestRays, s.closestHits);
        bind(s.layeredTable, s.layeredRays, s.layeredHits);
        bind(s.discoveryTable, s.seeds, s.vertices);
        dispatch(closest_, s.closestTable, s.batch.closestCount);
        dispatch(layered_, s.layeredTable, s.batch.layeredCount);
        dispatch(discovery_, s.discoveryTable, s.batch.discoveryCount);
        [encoder endEncoding];
        [s.commands endCommandBuffer];
        return true;
    }

    TraceCapacity capacity_;
    bool ready_ = false;
    id<MTLDevice> device_;
    id<MTL4CommandQueue> queue_;
    id<MTLSharedEvent> event_;
    id<MTLResidencySet> residency_;
    id<MTLComputePipelineState> closest_, layered_, discovery_;
    id<MTLBuffer> triangleFacets_, facets_;
    std::array<BatchStorage, kMaxBatchesInFlight> storage_;
    std::deque<BatchStorage*> free_, done_;
    std::deque<Retired> retired_;
    std::unordered_map<std::uint32_t, MeshEntry> entries_;
    std::vector<std::uint32_t> firstTriangle_, firstFacet_;
    const Scene* scene_ = nullptr;
    std::uint64_t meshRevision_ = 0, materialRevision_ = 0, submitted_ = 0;
};

}  // namespace

std::unique_ptr<TraceBackend> makeMetalTraceBackend(const MetalTraceConfig& config, std::string* error) {
    @autoreleasepool {
        auto backend = std::make_unique<MetalTraceBackend>(config, error);
        if (!backend->ready()) return nullptr;
        return backend;
    }
}

}  // namespace mspatial::acoustics
