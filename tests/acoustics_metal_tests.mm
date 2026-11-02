// GPU check of the Metal acoustic backend against the CPU backend. Needs the GPU, so it stays out of CTest: run it under the GPU lock.
//   marathon_spatial_acoustics_metal_tests <path to acoustic_trace.metallib>
// Touches no audio device.
#import <Foundation/Foundation.h>

#include "acoustics/AcousticSimulator.hpp"
#include "acoustics/MetalTraceBackend.hpp"

#include <chrono>
#include <cstdio>

using namespace mspatial::acoustics;

namespace {

int failures = 0;

void quad(std::vector<float>& positions, std::vector<std::uint32_t>& indices, std::vector<std::uint32_t>& materials, Vec3 a, Vec3 b, Vec3 c, Vec3 d, std::uint32_t material) {
    const std::uint32_t base = std::uint32_t(positions.size() / 3);
    for (Vec3 p : {a, b, c, d}) positions.insert(positions.end(), {p.x, p.y, p.z});
    indices.insert(indices.end(), {base, base + 1, base + 2, base, base + 2, base + 3});
    materials.insert(materials.end(), {material, material});
}

// Shoebox with a mid-room scaled and rotated partition, so instancing, normals and facets all take part.
std::uint32_t buildScene(Scene& scene) {
    Material wall;
    wall.absorption = uniformBands(0.15f);
    wall.scattering = uniformBands(0.3f);
    Material partition;
    partition.absorption = uniformBands(0.4f);
    partition.scattering = uniformBands(0.f);
    partition.transmission = uniformBands(0.05f);
    const std::uint32_t wallId = scene.addMaterial(wall), partitionId = scene.addMaterial(partition);
    std::vector<float> p;
    std::vector<std::uint32_t> i, m;
    const Vec3 lo = vec3(0, 0, 0), hi = vec3(10, 7, 3);
    const Vec3 c[8] = {vec3(lo.x, lo.y, lo.z), vec3(hi.x, lo.y, lo.z), vec3(hi.x, hi.y, lo.z), vec3(lo.x, hi.y, lo.z),
                       vec3(lo.x, lo.y, hi.z), vec3(hi.x, lo.y, hi.z), vec3(hi.x, hi.y, hi.z), vec3(lo.x, hi.y, hi.z)};
    quad(p, i, m, c[0], c[4], c[7], c[3], wallId);
    quad(p, i, m, c[1], c[2], c[6], c[5], wallId);
    quad(p, i, m, c[0], c[1], c[5], c[4], wallId);
    quad(p, i, m, c[3], c[7], c[6], c[2], wallId);
    quad(p, i, m, c[0], c[3], c[2], c[1], wallId);
    quad(p, i, m, c[4], c[5], c[6], c[7], wallId);
    scene.addInstance(scene.addMesh({p, i, m}), Transform{});
    std::vector<float> p2;
    std::vector<std::uint32_t> i2, m2;
    quad(p2, i2, m2, vec3(0, -1, 0), vec3(0, 1, 0), vec3(0, 1, 1), vec3(0, -1, 1), partitionId);
    Transform t;
    const float cosine = std::cos(0.5f), sine = std::sin(0.5f);
    t.rows[0][0] = cosine;
    t.rows[0][1] = -sine * 2.f;
    t.rows[0][3] = 5.f;
    t.rows[1][0] = sine;
    t.rows[1][1] = cosine * 2.f;
    t.rows[1][3] = 3.5f;
    t.rows[2][2] = 2.f;
    return scene.addInstance(scene.addMesh({p2, i2, m2}), t);
}

// Same rays through both backends.
void compare(TraceBackend& cpu, TraceBackend& gpu, const Scene& scene) {
    cpu.syncScene(scene);
    gpu.syncScene(scene);
    TraceBatch* a = cpu.acquireBatch();
    TraceBatch* b = gpu.acquireBatch();
    const std::uint32_t count = 1024;
    for (TraceBatch* batch : {a, b}) {
        batch->closestCount = batch->layeredCount = count / 4;
        batch->discoveryCount = 128;
    }
    for (std::uint32_t k = 0; k < count / 4; ++k) {
        const float z = 1.f - 2.f * (float(k) + 0.5f) / float(count / 4), r = std::sqrt(1.f - z * z), phi = float(k) * 2.399963f;
        const AcousticRayRecord ray{3.f, 2.f, 1.5f, 1e-3f, r * std::cos(phi), r * std::sin(phi), z, 100.f};
        a->closestRays[k] = b->closestRays[k] = a->layeredRays[k] = b->layeredRays[k] = ray;
    }
    for (std::uint32_t k = 0; k < 128; ++k) {
        const float z = 1.f - 2.f * (float(k) + 0.5f) / 128.f, r = std::sqrt(1.f - z * z), phi = float(k) * 2.399963f;
        const AcousticDiscoverySeed seed{3.f, 2.f, 1.5f, k * 2654435761u, r * std::cos(phi), r * std::sin(phi), z, 0};
        a->discoverySeeds[k] = b->discoverySeeds[k] = seed;
    }
    cpu.submit(a);
    gpu.submit(b);
    TraceBatch* doneA = nullptr;
    TraceBatch* doneB = nullptr;
    while (!(doneA = cpu.pollCompleted())) {}
    const auto start = std::chrono::steady_clock::now();
    while (!(doneB = gpu.pollCompleted())) {
        if (std::chrono::steady_clock::now() - start > std::chrono::seconds(10)) {
            std::printf("FAIL GPU batch did not finish in 10 s\n");
            ++failures;
            return;
        }
    }
    if (doneB->failed) {
        std::printf("FAIL GPU batch reported an error\n");
        ++failures;
        return;
    }
    std::uint32_t closestDiff = 0, layerDiff = 0, walkDiff = 0;
    for (std::uint32_t k = 0; k < count / 4; ++k) {
        const AcousticHitRecord& x = doneA->closestHits[k];
        const AcousticHitRecord& y = doneB->closestHits[k];
        closestDiff += x.instance != y.instance || x.facet != y.facet || (x.instance != ACOUSTIC_MISS && std::fabs(x.t - y.t) > 1e-3f);
        layerDiff += doneA->layeredHits[k].count != doneB->layeredHits[k].count;
    }
    const std::uint32_t maxOrder = cpu.capacity().maxOrder;
    for (std::uint32_t k = 0; k < 128; ++k) {
        for (std::uint32_t v = 0; v < maxOrder; ++v) {
            const AcousticDiscoveryVertex& x = doneA->discoveryVertices[k * maxOrder + v];
            const AcousticDiscoveryVertex& y = doneB->discoveryVertices[k * maxOrder + v];
            walkDiff += x.instance != y.instance || x.facet != y.facet;
            if (x.instance == ACOUSTIC_MISS) break;
        }
    }
    std::printf("%s closest differences %u, layered count differences %u, walk differences %u\n", closestDiff + layerDiff + walkDiff ? "FAIL" : "  ok", closestDiff, layerDiff, walkDiff);
    failures += closestDiff + layerDiff + walkDiff != 0;
    cpu.release(doneA);
    gpu.release(doneB);
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::printf("usage: %s <acoustic_trace.metallib>\n", argv[0]);
        return 2;
    }
    Scene scene;
    const std::uint32_t partition = buildScene(scene);
    MetalTraceConfig config;
    config.libraryPath = argv[1];
    std::string error;
    auto gpu = makeMetalTraceBackend(config, &error);
    if (!gpu) {
        std::printf("FAIL no Metal backend: %s\n", error.c_str());
        return 1;
    }
    auto cpu = makeCpuTraceBackend();
    compare(*cpu, *gpu, scene);
    // Move an instance: the next batch must use the new position.
    Transform moved = scene.instances()[partition].objectToWorld;
    moved.rows[1][3] -= 1.f;
    scene.setInstanceTransform(partition, moved);
    compare(*cpu, *gpu, scene);

    // The whole simulator on the GPU, ticking at 60 Hz without ever waiting for it.
    SimulatorConfig simulation;
    Simulator simulator(scene, *gpu, simulation);
    simulator.setListener(vec3(2, 2, 1.5f));
    for (int s = 0; s < 32; ++s) simulator.addSource({vec3(1.f + float(s % 8) * 1.1f, 1.f + float(s / 8) * 1.2f, 1.2f), 0.1f});
    double worst = 0.0, total = 0.0;
    for (int tick = 0; tick < 600; ++tick) {
        const auto begin = std::chrono::steady_clock::now();
        simulator.tick(1.f / 60.f);
        const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count();
        worst = std::max(worst, seconds);
        total += seconds;
        [NSThread sleepForTimeInterval:1.0 / 60.0];
    }
    std::uint32_t paths = 0;
    const AcousticSnapshot& snapshot = simulator.snapshots().read();
    for (std::uint32_t s = 0; s < 32; ++s)
        for (const ReflectionPath& path : snapshot.sources[s].reflections) paths += path.id != 0;
    std::printf("  measured: host tick %.3f ms mean, %.3f ms worst; %u paths over 32 sources\n", total / 600.0 * 1000.0, worst * 1000.0, paths);
    if (paths < 32 * 6) {
        std::printf("FAIL too few paths\n");
        ++failures;
    }
    std::printf(failures ? "FAILED\n" : "all passed\n");
    return failures ? 1 : 0;
}
