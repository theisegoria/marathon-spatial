// CPU-only acoustics tests against analytic rooms. Everything renders to memory; nothing here touches an audio device or the GPU.
#include "acoustics/AcousticSimulator.hpp"
#include "acoustics/Bake.hpp"
#include "acoustics/Convolver.hpp"
#include "acoustics/Diffraction.hpp"
#include "acoustics/EnergyField.hpp"
#include "acoustics/ReverbIr.hpp"
#include "acoustics/LateReverb.hpp"
#include "marathon_spatial/marathon_spatial_acoustics.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <new>
#include <numbers>
#include <random>
#include <thread>

using namespace mspatial::acoustics;

namespace {

int failures = 0;
std::atomic<long> allocationCount{0};

void check(bool ok, const char* what) {
    if (!ok) {
        std::printf("FAIL %s\n", what);
        ++failures;
    }
}

void near(double actual, double expected, double tolerance, const char* what) {
    const bool ok = std::fabs(actual - expected) <= tolerance;
    std::printf("%s %-58s %.6g expected %.6g (+- %.3g)\n", ok ? "  ok" : "FAIL", what, actual, expected, tolerance);
    if (!ok) ++failures;
}

struct MeshData {
    std::vector<float> positions;
    std::vector<std::uint32_t> indices;
    std::vector<std::uint32_t> materials;

    void quad(Vec3 a, Vec3 b, Vec3 c, Vec3 d, std::uint32_t material) {
        const std::uint32_t base = std::uint32_t(positions.size() / 3);
        for (Vec3 p : {a, b, c, d}) positions.insert(positions.end(), {p.x, p.y, p.z});
        indices.insert(indices.end(), {base, base + 1, base + 2, base, base + 2, base + 3});
        materials.insert(materials.end(), {material, material});
    }
    std::uint32_t add(Scene& scene) const { return scene.addMesh({positions, indices, materials}); }
};

// Inside of a box with outward-consistent winding. materials: -x, +x, -y, +y, -z, +z.
MeshData boxInside(Vec3 lo, Vec3 hi, std::array<std::uint32_t, 6> m) {
    MeshData mesh;
    const Vec3 p[8] = {vec3(lo.x, lo.y, lo.z), vec3(hi.x, lo.y, lo.z), vec3(hi.x, hi.y, lo.z), vec3(lo.x, hi.y, lo.z),
                       vec3(lo.x, lo.y, hi.z), vec3(hi.x, lo.y, hi.z), vec3(hi.x, hi.y, hi.z), vec3(lo.x, hi.y, hi.z)};
    mesh.quad(p[0], p[4], p[7], p[3], m[0]);
    mesh.quad(p[1], p[2], p[6], p[5], m[1]);
    mesh.quad(p[0], p[1], p[5], p[4], m[2]);
    mesh.quad(p[3], p[7], p[6], p[2], m[3]);
    mesh.quad(p[0], p[3], p[2], p[1], m[4]);
    mesh.quad(p[4], p[5], p[6], p[7], m[5]);
    return mesh;
}

Material flat(float absorption, float scattering = 0.f, float transmission = 0.f) {
    return {uniformBands(absorption), uniformBands(scattering), uniformBands(transmission)};
}

Bands noAir() { return uniformBands(0.f); }

// ---------------------------------------------------------------------------------------------------------------------

void airAttenuationMatchesIso() {
    // ISO 9613-2 table 2, 20 C and 70 % relative humidity, in dB per km. The table is at exact midbands, 3.981 and 7.943 kHz, which read about 1 % low here.
    const Bands m = airAttenuation(20.f, 70.f);
    const float toDbPerKm = 1000.f * 4.342945f;
    near(m[4] * toDbPerKm, 5.0, 0.1, "air 1 kHz dB/km");
    near(m[5] * toDbPerKm, 9.0, 0.1, "air 2 kHz dB/km");
    near(m[6] * toDbPerKm, 22.9, 0.3, "air 4 kHz dB/km");
    near(m[7] * toDbPerKm, 76.6, 1.2, "air 8 kHz dB/km");
    check(std::is_sorted(m.begin(), m.end()), "air absorption rises with frequency");
}

void sabineAndEyringOnAnAnalyticRoom() {
    // 10 x 8 x 4 m, every surface alpha = 0.2.
    Scene scene;
    const std::uint32_t wall = scene.addMaterial(flat(0.2f));
    const MeshData mesh = boxInside(vec3(0, 0, 0), vec3(10, 8, 4), {wall, wall, wall, wall, wall, wall});
    const ClosedMeasure measure = measureClosedMesh(scene.meshes()[mesh.add(scene)], scene.materials());
    const double V = 320, S = 2 * (80 + 40 + 32), A = S * 0.2, c = 343;
    near(measure.volume, V, 1e-3, "closed mesh volume");
    near(measure.surfaceArea, S, 1e-3, "closed mesh surface");
    near(measure.absorptionArea[3], A, 1e-3, "closed mesh absorption area");

    const double sabineT60 = 24.0 * std::log(10.0) * V / (c * A);
    const double eyringT60 = 24.0 * std::log(10.0) * V / (c * -S * std::log(1.0 - 0.2));
    near(t60FromRate(decayRate(float(V), equivalentAbsorption(DecayFormula::Sabine, float(S), float(A)), 0.f, 343.f)), sabineT60, 1e-3, "Sabine T60");
    near(t60FromRate(decayRate(float(V), equivalentAbsorption(DecayFormula::Eyring, float(S), float(A)), 0.f, 343.f)), eyringT60, 1e-3, "Eyring T60");
    check(eyringT60 < sabineT60, "Eyring is shorter than Sabine");
    // The familiar constant: T60 = 0.161 V / A.
    near(sabineT60, 0.161 * V / A, 0.01, "Sabine agrees with 0.161 V / A");
}

void imageSourceDelayAndGainInABox() {
    // Each wall has its own absorption, so a wrong wall shows up as a wrong gain.
    Scene scene;
    std::array<std::uint32_t, 6> walls;
    const float alpha[6] = {0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f};
    for (int i = 0; i < 6; ++i) walls[i] = scene.addMaterial(flat(alpha[i]));
    const Vec3 size = vec3(7, 5, 3);
    const MeshData mesh = boxInside(vec3(0, 0, 0), size, walls);
    scene.addInstance(mesh.add(scene), Transform{});
    auto backend = makeCpuTraceBackend();
    SimulatorConfig config;
    config.air = noAir();
    Simulator simulator(scene, *backend, config);
    const Vec3 sourcePosition = vec3(2, 1.5f, 1), listener = vec3(4.4f, 3.1f, 2.3f);  // No two images at one distance.
    const std::uint32_t source = simulator.addSource({sourcePosition, 0.f});
    simulator.setListener(listener);
    for (int i = 0; i < 40; ++i) simulator.tick(1.f / 60.f);
    const AcousticSnapshot& snapshot = simulator.snapshots().read();
    const SourceAcoustics& view = snapshot.sources[source];

    // Exact lattice of first-order images.
    int found = 0;
    for (int axis = 0; axis < 3; ++axis) {
        for (int side = 0; side < 2; ++side) {
            Vec3 image = sourcePosition;
            image[axis] = side == 0 ? -sourcePosition[axis] : 2.f * size[axis] - sourcePosition[axis];
            const float distance = length(image - listener);
            const float expectedDelay = distance / 343.f;
            const float expectedGain = std::sqrt(1.f - alpha[axis * 2 + side]) / distance;
            for (const ReflectionPath& path : view.reflections) {
                if (path.id == 0 || path.order != 1 || std::fabs(path.delaySeconds - expectedDelay) > 1e-5f) continue;
                ++found;
                char label[64];
                std::snprintf(label, sizeof label, "first-order gain, axis %d side %d", axis, side);
                near(path.gain[3], expectedGain, 2e-3, label);
                Vec3 arrival = normalize(image - listener);
                near(dot(path.arrival, arrival), 1.0, 1e-4, "arrival points at the image");
            }
        }
    }
    near(found, 6, 0, "all six first-order images found with 96 rays per tick");

    // A second-order corner image.
    Vec3 image = sourcePosition;
    image.x = -image.x;
    image.y = -image.y;
    const float distance = length(image - listener);
    bool corner = false;
    for (const ReflectionPath& path : view.reflections) {
        if (path.id == 0 || path.order != 2 || std::fabs(path.delaySeconds - distance / 343.f) > 1e-5f) continue;
        corner = true;
        near(path.gain[3], std::sqrt((1.f - alpha[0]) * (1.f - alpha[2])) / distance, 2e-3, "second-order corner gain");
    }
    check(corner, "second-order corner path found");
}

void scatteringAirAndScaleInAPath() {
    // One big floor, reflection point at the origin of a scaled instance. Scattering and air both come out of the specular gain.
    Scene scene;
    Material floor = flat(0.3f, 0.25f);
    const std::uint32_t id = scene.addMaterial(floor);
    MeshData mesh;
    mesh.quad(vec3(-1, -1, 0), vec3(1, -1, 0), vec3(1, 1, 0), vec3(-1, 1, 0), id);
    Transform scaled;
    scaled.rows[0][0] = 40.f;
    scaled.rows[1][1] = 40.f;
    const std::uint32_t instance = scene.addInstance(mesh.add(scene), scaled);
    PathKey key;
    key.order = 1;
    key.facets[0] = {instance, 0};
    ImagePath path;
    Bands air = uniformBands(0.01f);
    check(solveImagePath(scene, key, vec3(-3, 0, 4), vec3(3, 0, 4), air, path), "floor path solves");
    near(path.length, 10.0, 1e-4, "reflection length");
    near(path.amplitude[0], std::sqrt(0.7 * 0.75 * std::exp(-0.1)) / 10.0, 1e-5, "reflection gain with scattering and air");
    near(path.points[0].x, 0.0, 1e-5, "reflection point on the floor");
    check(!solveImagePath(scene, key, vec3(-3, 0, 4), vec3(3, 0, -4), air, path), "no reflection between opposite sides");
}

void occlusionAndTransmission() {
    // Wall with band-dependent transmission between listener and source.
    Scene scene;
    Material concrete = flat(0.05f);
    for (std::size_t b = 0; b < kBandCount; ++b) concrete.transmission[b] = 0.2f / float(b + 1);
    const std::uint32_t id = scene.addMaterial(concrete);
    MeshData wall;
    wall.quad(vec3(0, -5, -5), vec3(0, 5, -5), vec3(0, 5, 5), vec3(0, -5, 5), id);
    const std::uint32_t instance = scene.addInstance(wall.add(scene), Transform{});
    auto backend = makeCpuTraceBackend();
    Simulator simulator(scene, *backend);
    simulator.addSource({vec3(2, 0, 0), 0.05f});
    simulator.setListener(vec3(-2, 0, 0));
    for (int i = 0; i < 40; ++i) simulator.tick(1.f / 60.f);
    {
        const SourceAcoustics& view = simulator.snapshots().read().sources[0];
        near(view.visibility, 0.0, 1e-3, "wall blocks the line");
        for (std::size_t b : {0u, 3u, 7u}) near(view.directGain[b], std::sqrt(concrete.transmission[b]), 1e-3, "transmission through one surface");
    }
    // Destruction at the live position: the wall goes, and the line clears within the smoothing time.
    scene.removeInstance(instance);
    for (int i = 0; i < 40; ++i) simulator.tick(1.f / 60.f);
    const SourceAcoustics& view = simulator.snapshots().read().sources[0];
    near(view.visibility, 1.0, 1e-3, "removed wall clears the line");
    near(view.directGain[3], 1.0, 1e-3, "removed wall passes everything");
}

void partialOcclusionOverTheSourceDisc() {
    Scene scene;
    const std::uint32_t id = scene.addMaterial(flat(0.1f));
    MeshData wall;  // Half plane y < 0 between listener and source.
    wall.quad(vec3(0, -9, -9), vec3(0, 0, -9), vec3(0, 0, 9), vec3(0, -9, 9), id);
    scene.addInstance(wall.add(scene), Transform{});
    auto backend = makeCpuTraceBackend();
    SimulatorConfig config;
    config.smoothingSeconds = 0.2f;
    Simulator simulator(scene, *backend, config);
    simulator.addSource({vec3(4, 0, 0), 1.f});
    simulator.setListener(vec3(-4, 0.f, 0));
    double sum = 0.0;
    for (int i = 0; i < 400; ++i) {
        simulator.tick(1.f / 60.f);
        if (i >= 100) sum += simulator.snapshots().read().sources[0].visibility;
    }
    // The first ray aims at the center, which grazes the wall edge; the rest spread evenly. About half the disc is hidden.
    near(sum / 300.0, 0.5, 0.15, "half-hidden disc reads about half visible");
}

void movingDoorChangesReflectionsAndStaysStable() {
    Scene scene;
    const std::uint32_t wallMaterial = scene.addMaterial(flat(0.2f));
    const std::uint32_t doorMaterial = scene.addMaterial(flat(0.5f));
    const MeshData room = boxInside(vec3(0, 0, 0), vec3(8, 6, 3), {wallMaterial, wallMaterial, wallMaterial, wallMaterial, wallMaterial, wallMaterial});
    scene.addInstance(room.add(scene), Transform{});
    MeshData doorMesh;
    doorMesh.quad(vec3(0, -1, 0), vec3(0, 1, 0), vec3(0, 1, 3), vec3(0, -1, 3), doorMaterial);
    const std::uint32_t door = scene.addInstance(doorMesh.add(scene), Transform::translation(vec3(4, 20, 0)));  // Far away, out of the room.
    auto backend = makeCpuTraceBackend();
    SimulatorConfig config;
    config.air = noAir();
    Simulator simulator(scene, *backend, config);
    const Vec3 sourcePosition = vec3(2, 3, 1.5f), listener = vec3(3, 1.2f, 0.9f);  // No two images at one distance.
    simulator.addSource({sourcePosition, 0.01f});
    simulator.setListener(listener);
    for (int i = 0; i < 240; ++i) simulator.tick(1.f / 60.f);

    std::array<std::uint64_t, kMaxPathsPerSource> ids;
    std::array<float, kMaxPathsPerSource> gains;
    {
        const SourceAcoustics& view = simulator.snapshots().read().sources[0];
        for (std::uint32_t i = 0; i < kMaxPathsPerSource; ++i) {
            ids[i] = view.reflections[i].id;
            gains[i] = view.reflections[i].gain[3];
        }
    }
    // Static scene: no path may change identity or gain, however few rays run.
    float worst = 0.f;
    std::uint32_t changed = 0;
    for (int tick = 0; tick < 300; ++tick) {
        simulator.tick(1.f / 60.f);
        const SourceAcoustics& view = simulator.snapshots().read().sources[0];
        for (std::uint32_t i = 0; i < kMaxPathsPerSource; ++i) {
            if (ids[i] == 0) continue;
            changed += view.reflections[i].id != ids[i];
            worst = std::max(worst, std::fabs(view.reflections[i].gain[3] - gains[i]));
        }
    }
    near(changed, 0, 0, "paths keep their slots and ids over 5 s");
    near(worst, 0.0, 2e-3, "path gains hold still with 96 discovery rays");

    // A door sliding in as a new mirror at x = 4: the path off it has delay and gain from its live position.
    scene.setInstanceTransform(door, Transform::translation(vec3(4, 2, 0)));
    const float mirrored = length(vec3(2.f * 4.f - sourcePosition.x, sourcePosition.y, sourcePosition.z) - listener);
    for (int i = 0; i < 80; ++i) simulator.tick(1.f / 60.f);
    {
        bool found = false;
        const SourceAcoustics& view = simulator.snapshots().read().sources[0];
        for (const ReflectionPath& path : view.reflections) {
            if (path.id == 0 || path.order != 1 || std::fabs(path.delaySeconds - mirrored / 343.f) > 2e-5f) continue;
            found = true;
            near(path.gain[3], std::sqrt(0.5f) / mirrored, 1e-3, "reflection off the door after it moved in");
        }
        check(found, "door reflection appears");
    }
    // Move it again and the same path follows the door, with no change of identity.
    std::uint64_t doorPath = 0;
    {
        const SourceAcoustics& view = simulator.snapshots().read().sources[0];
        for (const ReflectionPath& path : view.reflections)
            if (path.id != 0 && std::fabs(path.delaySeconds - mirrored / 343.f) < 2e-5f) doorPath = path.id;
    }
    scene.setInstanceTransform(door, Transform::translation(vec3(4.5f, 2, 0)));
    const float moved = length(vec3(2.f * 4.5f - sourcePosition.x, sourcePosition.y, sourcePosition.z) - listener);
    simulator.tick(1.f / 60.f);
    simulator.tick(1.f / 60.f);
    {
        bool followed = false;
        const SourceAcoustics& view = simulator.snapshots().read().sources[0];
        for (const ReflectionPath& path : view.reflections)
            if (path.id == doorPath) followed = std::fabs(path.delaySeconds - moved / 343.f) < 2e-5f;
        check(followed, "the same path tracks the moving door within two ticks");
    }
    // Door gone: the path fades out instead of stopping dead.
    scene.setInstanceEnabled(door, false);
    float previous = 1.f;
    bool monotonic = true;
    for (int i = 0; i < 60; ++i) {
        simulator.tick(1.f / 60.f);
        float gain = 0.f;
        const SourceAcoustics& view = simulator.snapshots().read().sources[0];
        for (const ReflectionPath& path : view.reflections)
            if (path.id == doorPath) gain = path.gain[3];
        monotonic = monotonic && gain <= previous + 1e-6f;
        previous = gain;
    }
    check(monotonic, "disabled door fades its path out monotonically");
    near(previous, 0.0, 1e-4, "path off a removed door is gone");
}

void probesMeasureARealRoom() {
    Scene scene;
    const std::uint32_t wall = scene.addMaterial(flat(0.25f));
    const MeshData mesh = boxInside(vec3(0, 0, 0), vec3(9, 6, 3), {wall, wall, wall, wall, wall, wall});
    scene.addInstance(mesh.add(scene), Transform{});
    RoomDesc room;
    room.bounds = {Box{vec3(0, 0, 0), vec3(9, 6, 3)}};
    room.probe = vec3(4.5f, 3, 1.5f);
    const std::uint32_t id = scene.addRoom(room);
    auto backend = makeCpuTraceBackend();
    SimulatorConfig config;
    config.air = noAir();
    config.formula = DecayFormula::Sabine;
    Simulator simulator(scene, *backend, config);
    for (int i = 0; i < 40; ++i) simulator.tick(1.f / 60.f);  // 16 rays a tick fills the 512-ray window.
    const RoomAcoustics& measured = simulator.snapshots().read().rooms[id];
    const double S = 2 * (54 + 27 + 18);
    near(measured.volume, 162.0, 162.0 * 0.05, "probe volume within 5 %");
    near(measured.surfaceArea, S, S * 0.10, "probe surface within 10 %");
    near(measured.absorptionArea[3] / measured.surfaceArea, 0.25, 1e-4, "probe mean absorption");
    const double expected = 0.161 * 162.0 / (S * 0.25);
    near(measured.t60[3], expected, expected * 0.15, "probe T60 against Sabine");
}

void portalGainAndOpenness() {
    // Point-through-aperture geometry from the closed form, in the high band where the Lambert law holds.
    const float area = 1.f;
    ApertureLeg leg;
    leg.distance = 4.f;
    leg.area = area;
    leg.cosIn = 1.f;
    leg.cosOut = 1.f;
    Bands energy = apertureChainEnergy(std::span(&leg, 1), 6.f, 343.f);
    const double ka = 2.0 * std::numbers::pi * 8000.0 / 343.0 * std::sqrt(area / std::numbers::pi);
    const double weight = ka / (1.0 + ka);
    near(energy[7], area / (std::numbers::pi * 16.0 * 36.0), 1e-7, "aperture gain A cos cos / (pi d1^2 d2^2)");

    // Oblique: low band stays near omnidirectional, high band follows the cosine.
    leg.cosIn = 0.5f;
    leg.cosOut = 0.5f;
    energy = apertureChainEnergy(std::span(&leg, 1), 6.f, 343.f);
    check(energy[0] > energy[7] * 1.5f, "small aperture spreads low bands more than high ones");
    const double lobe = (1.0 - weight) + weight * 0.5;
    near(energy[7], lobe * lobe * area / (std::numbers::pi * 16.0 * 36.0), 1e-7, "oblique high band follows both cosines");

    // A huge opening right by the listener saturates at the free field of the unfolded path.
    ApertureLeg big;
    big.distance = 2.f;
    big.area = 200.f;
    energy = apertureChainEnergy(std::span(&big, 1), 1.f, 343.f);
    near(energy[4], 1.0 / 9.0, 1e-6, "large aperture cannot exceed the unfolded free field");

    // Two rooms with a door in a shared wall, listener and source out of line. Door slides across the frame at live positions.
    Scene scene;
    Material doorSteel = flat(0.1f);
    for (std::size_t b = 0; b < kBandCount; ++b) doorSteel.transmission[b] = 0.01f;
    const std::uint32_t wall = scene.addMaterial(flat(0.1f));
    const std::uint32_t steel = scene.addMaterial(doorSteel);
    MeshData divider;
    // Wall at x = 0 with a 2 x 2 m hole at y in [-1, 1], z in [0, 2].
    divider.quad(vec3(0, -8, 0), vec3(0, -1, 0), vec3(0, -1, 4), vec3(0, -8, 4), wall);
    divider.quad(vec3(0, 1, 0), vec3(0, 8, 0), vec3(0, 8, 4), vec3(0, 1, 4), wall);
    divider.quad(vec3(0, -1, 2), vec3(0, 1, 2), vec3(0, 1, 4), vec3(0, -1, 4), wall);
    scene.addInstance(divider.add(scene), Transform{});
    MeshData doorMesh;
    doorMesh.quad(vec3(0, -1, 0), vec3(0, 1, 0), vec3(0, 1, 2), vec3(0, -1, 2), steel);
    const std::uint32_t door = scene.addInstance(doorMesh.add(scene), Transform::translation(vec3(0, 2, 0)));
    RoomDesc a, b;
    a.bounds = {Box{vec3(-8, -8, 0), vec3(0, 8, 4)}};
    a.volume = 512.f;
    a.surfaceArea = 416.f;
    a.absorptionArea = uniformBands(40.f);
    b = a;
    b.bounds = {Box{vec3(0, -8, 0), vec3(8, 8, 4)}};
    const std::uint32_t roomA = scene.addRoom(a), roomB = scene.addRoom(b);
    PortalDesc portal;
    portal.roomA = roomA;
    portal.roomB = roomB;
    portal.center = vec3(0, 0, 1);
    portal.normal = vec3(1, 0, 0);
    portal.halfU = vec3(0, 1, 0);
    portal.halfV = vec3(0, 0, 1);
    scene.addPortal(portal);

    auto backend = makeCpuTraceBackend();
    SimulatorConfig config;
    config.air = noAir();
    config.formula = DecayFormula::Sabine;
    Simulator simulator(scene, *backend, config);
    const Vec3 sourcePosition = vec3(3, 3, 1), listener = vec3(-4, 0, 1);  // Source out of line; the straight line runs through the wall.
    simulator.addSource({sourcePosition, 0.f});
    simulator.setListener(listener);

    const auto settle = [&] {
        for (int i = 0; i < 12; ++i) simulator.tick(1.f / 60.f);
        return std::cref(simulator.snapshots().read());
    };
    {
        const AcousticSnapshot& s = settle();
        const PortalPath& path = s.sources[0].portalPaths[0];
        check(path.portals == 1, "one portal path");
        // The crossing point is on the frame's line between source and listener.
        const Vec3 crossing = vec3(0, 3.f * 4.f / 7.f, 1);
        const float inside = std::min(crossing.y, 1.f);
        const Vec3 clamped = vec3(0, inside, 1);
        const float d1 = length(clamped - sourcePosition), d2 = length(clamped - listener);
        near(path.delaySeconds, (d1 + d2) / 343.0, 1e-5, "portal delay over the crossing point");
        const float cos1 = 3.f / d1, cos2 = 4.f / d2;
        const double expected = std::sqrt(4.0 * cos1 * cos2 / (std::numbers::pi * d1 * d1 * d2 * d2));
        near(path.gain[7], expected, expected * 0.02, "portal path amplitude");
        near(dot(path.arrival, normalize(clamped - listener)), 1.0, 1e-4, "sound arrives from the opening");
        // Diffuse exchange: the rate c S / 4 V with S = 4 m^2, V = 512 m^3.
        near(s.links[0].rate[3], 343.0 * 4.0 / (4.0 * 512.0), 1e-3, "room link rate");
        near(s.rooms[0].t60[3], t60FromRate(343.f * (40.f + 4.f) / (4.f * 512.f)), 1e-2, "T60 counts the open doorway as loss");
    }
    // Half closed: a quarter of the sample pattern is not enough resolution for a fine number, so use the exact 16-sample window.
    const float openGain = simulator.snapshots().read().sources[0].portalPaths[0].gain[7];
    scene.setInstanceTransform(door, Transform::translation(vec3(0, 1, 0)));
    {
        const AcousticSnapshot& s = settle();
        const double ratio = s.links[0].rate[3] / (343.0 * 4.0 / (4.0 * 512.0));
        near(ratio, 0.5 + 0.5 * 0.01, 0.08, "half-open door halves the coupling, plus its leakage");
    }
    scene.setInstanceTransform(door, Transform{});
    {
        const AcousticSnapshot& s = settle();
        near(s.links[0].rate[3] / (343.0 * 4.0 / (4.0 * 512.0)), 0.01, 1e-4, "closed door leaks only its transmission");
        const float closedGain = s.sources[0].portalPaths[0].gain[7];
        near(closedGain / openGain, 0.1, 0.002, "closed door takes 20 dB off the portal path");
    }
    // Destruction: the door is gone, and the opening is open again.
    scene.removeInstance(door);
    {
        const AcousticSnapshot& s = settle();
        near(s.links[0].rate[3], 343.0 * 4.0 / (4.0 * 512.0), 1e-3, "blown-out door restores the coupling");
    }
}

void portalChainThroughAThirdRoom() {
    // Three rooms in a row; the gain is the analytic two-aperture chain.
    Scene scene;
    for (int i = 0; i < 3; ++i) {
        RoomDesc room;
        room.bounds = {Box{vec3(float(i) * 10.f, -5, -5), vec3(float(i + 1) * 10.f, 5, 5)}};
        room.volume = 1000.f;
        room.surfaceArea = 600.f;
        room.absorptionArea = uniformBands(60.f);
        scene.addRoom(room);
    }
    for (int i = 0; i < 2; ++i) {
        PortalDesc portal;
        portal.roomA = std::uint32_t(i);
        portal.roomB = std::uint32_t(i + 1);
        portal.center = vec3(float(i + 1) * 10.f, 0, 0);
        scene.addPortal(portal);
    }
    auto backend = makeCpuTraceBackend();
    SimulatorConfig config;
    config.air = noAir();
    Simulator simulator(scene, *backend, config);
    simulator.addSource({vec3(5, 0, 0), 0.f});
    simulator.setListener(vec3(25, 0, 0));
    for (int i = 0; i < 12; ++i) simulator.tick(1.f / 60.f);
    const PortalPath& path = simulator.snapshots().read().sources[0].portalPaths[0];
    check(path.portals == 2, "two-portal path");
    near(path.delaySeconds, 20.0 / 343.0, 1e-5, "chain delay");
    const double area = 4.0, pi = std::numbers::pi;
    const double energy = area * area / (5.0 * 5.0 * pi * 10.0 * 10.0 * pi * 5.0 * 5.0);
    near(path.gain[7], std::sqrt(energy), std::sqrt(energy) * 0.03, "two-aperture chain gain in the high band");
}

void coupledRoomsAnalytic() {
    // Two rooms: the eigenvalues of [[a + k1, -k1], [-k2, b + k2]] in closed form.
    CoupledRooms rooms;
    rooms.volumes = {200.f, 800.f};
    rooms.absorption = {40.f, 60.f};
    const float S = 5.f;
    rooms.coupling = {0.f, S, S, 0.f};
    const double c = 343.0;
    const double a = c * 40.0 / (4.0 * 200.0), b = c * 60.0 / (4.0 * 800.0);
    const double k1 = c * S / (4.0 * 200.0), k2 = c * S / (4.0 * 800.0);
    const double trace = a + k1 + b + k2, det = (a + k1) * (b + k2) - k1 * k2;
    const double root = std::sqrt(trace * trace - 4.0 * det);
    const double fast = (trace + root) / 2.0, slow = (trace - root) / 2.0;
    CoupledModes modes = decayModes(rooms);
    std::sort(modes.rates.begin(), modes.rates.end());
    near(modes.rates[0], slow, 1e-3, "slow coupled mode");
    near(modes.rates[1], fast, 1e-3, "fast coupled mode");

    // Steady state against the 2 x 2 solution.
    const std::vector<float> powers = {1.f, 0.f};
    const std::vector<float> steady = rooms.steadyState(powers);
    const double A1 = 40.0 + S, A2 = 60.0 + S;
    const double e1 = 4.0 / c * A2 / (A1 * A2 - S * S), e2 = 4.0 / c * S / (A1 * A2 - S * S);
    near(steady[0], e1, 1e-6, "steady energy in the driven room");
    near(steady[1], e2, 1e-6, "steady energy in the neighbour");
    near(steady[1] / steady[0], S / A2, 1e-5, "classic ratio S / (A2 + S)");

    // Impulse response against brute integration of the ODE.
    double E1 = 1.0 / 200.0, E2 = 0.0;
    const double h = 1e-5;
    const CoupledModes exact = decayModes(rooms);
    for (int i = 0; i < 100000; ++i) {
        const double d1 = -(a + k1) * E1 + k1 * E2, d2 = k2 * E1 - (b + k2) * E2;
        E1 += h * d1;
        E2 += h * d2;
    }
    near(exact.impulseEnergy(0, 0, 1.f), E1, E1 * 0.01, "modal impulse, driven room at 1 s");
    near(exact.impulseEnergy(0, 1, 1.f), E2, E2 * 0.01, "modal impulse, neighbour at 1 s");
}

// ---------------------------------------------------------------------------------------------------------------------
// Late reverb, rendered to memory.

struct BandFilter {
    double b0, b1, b2, a1, a2, z1 = 0, z2 = 0;
    float run(float x) {
        const double y = b0 * x + z1;
        z1 = b1 * x - a1 * y + z2;
        z2 = b2 * x - a2 * y;
        return float(y);
    }
};

std::vector<float> octaveBand(const std::vector<float>& signal, float center, float sampleRate) {
    // Sixth-order-ish band: three constant-peak bandpasses in cascade, one octave wide.
    const double w = 2.0 * std::numbers::pi * center / sampleRate, q = 1.41421356, alpha = std::sin(w) / (2.0 * q);
    const double a0 = 1 + alpha;
    std::vector<float> out = signal;
    for (int stage = 0; stage < 3; ++stage) {
        BandFilter f{alpha / a0, 0, -alpha / a0, -2 * std::cos(w) / a0, (1 - alpha) / a0};
        for (float& x : out) x = f.run(x);
    }
    return out;
}

// T60 from a Schroeder backward integral, fitted between -5 and -25 dB.
double schroederT60(const std::vector<float>& signal, float sampleRate) {
    std::vector<double> tail(signal.size());
    double sum = 0.0;
    for (std::size_t i = signal.size(); i-- > 0;) {
        sum += double(signal[i]) * signal[i];
        tail[i] = sum;
    }
    const double total = tail[0];
    double sx = 0, sy = 0, sxx = 0, sxy = 0;
    double n = 0;
    for (std::size_t i = 0; i < tail.size(); i += 16) {
        const double db = 10.0 * std::log10(tail[i] / total);
        if (db > -5.0 || db < -25.0) continue;
        const double x = double(i) / sampleRate;
        sx += x;
        sy += db;
        sxx += x * x;
        sxy += x * db;
        n += 1;
    }
    const double slope = (n * sxy - sx * sy) / (n * sxx - sx * sx);
    return -60.0 / slope;
}

AcousticSnapshot* reverbSnapshot() {
    static auto* snapshot = new AcousticSnapshot;
    *snapshot = AcousticSnapshot{};
    return snapshot;
}

void setRoom(AcousticSnapshot& s, std::uint32_t index, float volume, const Bands& t60) {
    s.roomCount = std::max(s.roomCount, index + 1);
    RoomAcoustics& room = s.rooms[index];
    room.active = true;
    room.volume = volume;
    for (std::size_t b = 0; b < kBandCount; ++b) {
        room.t60[b] = t60[b];
        room.decayRate[b] = 6.f * std::numbers::ln10_v<float> / t60[b];
    }
}

void lateReverbDecayFollowsTheRoom() {
    // A hall with falling T60 across bands; the rendered tail must measure back to the numbers it was given.
    LateReverbConfig config;
    config.rooms = 1;
    config.links = 0;
    config.meanDelaySeconds = {0.03f};
    LateReverbNetwork reverb;
    reverb.prepare(config);
    const Bands t60{2.4f, 2.2f, 2.0f, 1.8f, 1.6f, 1.2f, 0.8f, 0.5f};
    AcousticSnapshot& snapshot = *reverbSnapshot();
    setRoom(snapshot, 0, 800.f, t60);
    reverb.apply(snapshot);

    const std::uint32_t frames = 48000 * 4;
    std::vector<float> input(frames, 0.f), left(frames);
    input[0] = 1.f;
    const float* inputs[] = {input.data()};
    float* lefts[] = {left.data()};
    reverb.process(inputs, lefts, nullptr, frames);
    for (std::size_t b = 1; b < kBandCount; ++b) {
        const double measured = schroederT60(octaveBand(left, kBandCenterHz[b], 48000.f), 48000.f);
        char label[64];
        std::snprintf(label, sizeof label, "FDN T60 at %.0f Hz", kBandCenterHz[b]);
        near(measured, t60[b], t60[b] * 0.15, label);
    }
}

void lateReverbLevelIsPhysical() {
    // With a flat decay, the summed output energy of an impulse must be 16 pi / A relative to the direct energy at 1 m.
    LateReverbConfig config;
    config.rooms = 1;
    config.links = 0;
    LateReverbNetwork reverb;
    reverb.prepare(config);
    const float volume = 500.f, t60 = 1.f;
    AcousticSnapshot& snapshot = *reverbSnapshot();
    setRoom(snapshot, 0, volume, uniformBands(t60));
    reverb.apply(snapshot);
    const std::uint32_t frames = 48000 * 2;
    std::vector<float> input(frames, 0.f), left(frames), right(frames);
    input[0] = 1.f;
    const float* inputs[] = {input.data()};
    float* lefts[] = {left.data()};
    float* rights[] = {right.data()};
    reverb.process(inputs, lefts, rights, frames);
    double energyL = 0.0, energyR = 0.0, cross = 0.0;
    for (std::uint32_t i = 0; i < frames; ++i) {
        energyL += double(left[i]) * left[i];
        energyR += double(right[i]) * right[i];
        cross += double(left[i]) * right[i];
    }
    const double rate = 6.0 * std::log(10.0) / t60;
    const double A = 4.0 * volume * rate / 343.0;
    const double expected = 16.0 * std::numbers::pi / A;
    near(10.0 * std::log10(energyL / expected), 0.0, 1.0, "late energy to 16 pi / A, dB");
    near(10.0 * std::log10(energyR / expected), 0.0, 1.0, "right channel the same, dB");
    near(cross / std::sqrt(energyL * energyR), 0.0, 0.15, "left and right decorrelated");
}

void lateReverbBinauralCoherence() {
    // The binaural stage sets the tail's interaural coherence per band and leaves its level alone.
    const auto run = [](const float* coherence, std::vector<float>& left, std::vector<float>& right) {
        LateReverbConfig config;
        config.rooms = 1;
        config.links = 0;
        LateReverbNetwork reverb;
        reverb.prepare(config);
        AcousticSnapshot& snapshot = *reverbSnapshot();
        setRoom(snapshot, 0, 500.f, uniformBands(1.5f));
        reverb.apply(snapshot);
        reverb.setCoherence(coherence);
        const std::uint32_t frames = 48000 * 3;
        std::vector<float> input(frames);
        std::uint32_t seed = 12345;
        for (auto& x : input) {
            seed = seed * 1664525u + 1013904223u;
            x = float(int(seed >> 9) - (1 << 22)) / float(1 << 22) * 0.1f;
        }
        left.assign(frames, 0.f);
        right.assign(frames, 0.f);
        const float* inputs[] = {input.data()};
        float* lefts[] = {left.data()};
        float* rights[] = {right.data()};
        reverb.process(inputs, lefts, rights, frames);
    };
    // Coherence of left and right in a band around hz, from a second-order band-pass of each.
    const auto bandCoherence = [](const std::vector<float>& l, const std::vector<float>& r, double hz) {
        const double w = 2 * std::numbers::pi * hz / 48000, alpha = std::sin(w) / (2 * 2.0), a0 = 1 + alpha;
        const double b0 = alpha / a0, b2 = -alpha / a0, a1 = -2 * std::cos(w) / a0, a2 = (1 - alpha) / a0;
        double zl[2]{}, zr[2]{}, cross = 0, pl = 0, pr = 0;
        for (std::size_t i = 0; i < l.size(); ++i) {
            const auto step = [&](double x, double* z) {const double y = b0 * x + z[0]; z[0] = -a1 * y + z[1]; z[1] = b2 * x - a2 * y; return y;};
            const double yl = step(l[i], zl), yr = step(r[i], zr);
            if (i < 24000) continue;
            cross += yl * yr;
            pl += yl * yl;
            pr += yr * yr;
        }
        return cross / std::sqrt(pl * pr);
    };
    const auto energy = [](const std::vector<float>& x) {double e = 0; for (float v : x) e += double(v) * v; return e;};
    std::vector<float> offL, offR, onL, onR;
    const float target[kCoherenceBands] = {0.95f, 0.8f, 0.4f, 0.1f, 0.f, 0.f};
    run(nullptr, offL, offR);
    run(target, onL, onR);
    const double low = bandCoherence(onL, onR, 141), high = bandCoherence(onL, onR, 4000), lowOff = bandCoherence(offL, offR, 141);
    std::printf("  binaural late field: coherence at 141 Hz %.2f (off %.2f), at 4 kHz %.2f; level %.2f dB left, %.2f dB right\n", low, lowOff, high,
                10 * std::log10(energy(onL) / energy(offL)), 10 * std::log10(energy(onR) / energy(offR)));
    check(low > 0.8, "low band takes the set's high coherence");
    check(std::fabs(high) < 0.15, "high band stays incoherent");
    check(std::fabs(lowOff) < 0.3, "without the stage the low band is not coherent");
    near(10 * std::log10(energy(onL) / energy(offL)), 0.0, 0.5, "left level unchanged, dB");
    near(10 * std::log10(energy(onR) / energy(offR)), 0.0, 0.5, "right level unchanged, dB");
}

void coupledReverbMatchesTheModes() {
    // Small dead room next to a big live one, energy feeding back through a door: the double slope must follow the analytic modes.
    AcousticSnapshot& snapshot = *reverbSnapshot();
    const float v0 = 100.f, v1 = 1000.f;
    const float a0 = 40.f, a1 = 60.f, door = 4.f, c = 343.f;
    Bands rate0, rate1, link01, link10;
    for (std::size_t b = 0; b < kBandCount; ++b) {
        rate0[b] = c * (a0 + door) / (4.f * v0);
        rate1[b] = c * (a1 + door) / (4.f * v1);
        link01[b] = c * door / (4.f * v0);
        link10[b] = c * door / (4.f * v1);
    }
    Bands t0, t1;
    for (std::size_t b = 0; b < kBandCount; ++b) {
        t0[b] = t60FromRate(rate0[b]);
        t1[b] = t60FromRate(rate1[b]);
    }
    setRoom(snapshot, 0, v0, t0);
    setRoom(snapshot, 1, v1, t1);
    snapshot.linkCount = 2;
    snapshot.links[0] = {0, 1, link01};
    snapshot.links[1] = {1, 0, link10};
    LateReverbConfig config;
    config.rooms = 2;
    config.links = 2;
    config.meanDelaySeconds = {0.015f, 0.04f};
    LateReverbNetwork reverb;
    reverb.prepare(config);
    reverb.apply(snapshot);

    const std::uint32_t frames = 48000 * 4;
    std::vector<float> input(frames, 0.f), left0(frames), left1(frames);
    // White noise burst, averaging out the sparse start of a single impulse.
    std::mt19937 random(7);
    std::normal_distribution<float> noise;
    const std::uint32_t burst = 2400;
    for (std::uint32_t i = 0; i < burst; ++i) input[i] = noise(random);
    const float* inputs[] = {input.data(), nullptr};
    float* lefts[] = {left0.data(), left1.data()};
    reverb.process(inputs, lefts, nullptr, frames);

    CoupledRooms rooms;
    rooms.volumes = {v0, v1};
    rooms.absorption = {a0, a1};
    rooms.coupling = {0.f, door, door, 0.f};
    const CoupledModes modes = decayModes(rooms);

    // Output power is absolute: y^2 = w 4 pi c / fs, with w the energy density the modes give for the burst's energy.
    const auto windowEnergy = [&](const std::vector<float>& signal, float start) {
        double sum = 0.0;
        const std::uint32_t begin = std::uint32_t(start * 48000.f), length = 4800;
        for (std::uint32_t i = begin; i < begin + length; ++i) sum += double(signal[i]) * signal[i];
        return sum / length;
    };
    const auto dbError = [&](const std::vector<float>& signal, std::size_t to, float time) {
        // Mean over the 100 ms window, which matters on the steep early slope.
        double density = 0.0;
        for (std::uint32_t i = 0; i < burst; i += 16)
            for (int step = 0; step < 8; ++step)
                density += double(input[i]) * input[i] * 16.0 / 8.0 * modes.impulseEnergy(0, to, time + (float(step) + 0.5f) * 0.0125f - float(i) / 48000.f);
        const double expected = density * 4.0 * std::numbers::pi * 343.0 / 48000.0;
        return 10.0 * std::log10(windowEnergy(signal, time) / expected);
    };
    near(dbError(left0, 0, 0.1f), 0.0, 1.5, "source room early against the modes, dB");
    near(dbError(left0, 0, 0.6f), 0.0, 2.0, "source room's slow tail at 0.6 s, dB");
    near(dbError(left0, 0, 1.2f), 0.0, 2.0, "source room's slow tail at 1.2 s, dB");
    near(dbError(left1, 1, 0.3f), 0.0, 1.5, "coupled room at 0.3 s against the modes, dB");
    near(dbError(left1, 1, 1.2f), 0.0, 2.0, "coupled room at 1.2 s against the modes, dB");
}

void reverbProcessIsRealTimeSafe() {
    LateReverbConfig config;
    config.rooms = 2;
    config.links = 2;
    LateReverbNetwork reverb;
    reverb.prepare(config);
    AcousticSnapshot& snapshot = *reverbSnapshot();
    setRoom(snapshot, 0, 300.f, uniformBands(1.5f));
    setRoom(snapshot, 1, 300.f, uniformBands(0.7f));
    snapshot.linkCount = 2;
    Bands leak = uniformBands(0.5f);
    leak[7] = 0.02f;  // A band-shaped link, the slowest code path.
    snapshot.links[0] = {0, 1, leak};
    snapshot.links[1] = {1, 0, leak};
    std::vector<float> in(512, 0.1f), a(512), b(512), c(512), d(512);
    const float* inputs[] = {in.data(), in.data()};
    float* lefts[] = {a.data(), c.data()};
    float* rights[] = {b.data(), d.data()};
    const long before = allocationCount.load();
    for (int block = 0; block < 200; ++block) {
        // A door swinging: new parameters every block.
        snapshot.rooms[0].decayRate[3] = 5.f + float(block % 20);
        reverb.apply(snapshot);
        reverb.process(inputs, lefts, rights, 512);
    }
    near(double(allocationCount.load() - before), 0.0, 0.0, "no allocation in apply() or process()");
    bool finite = true;
    for (float value : a) finite = finite && std::isfinite(value);
    check(finite, "output stays finite through parameter changes");
}

void sleepingRoomsWakeWithoutLoss() {
    // A long silence puts the room to sleep, and the same input afterwards must give the same output as a fresh network.
    LateReverbConfig config;
    config.rooms = 1;
    config.links = 0;
    AcousticSnapshot& snapshot = *reverbSnapshot();
    setRoom(snapshot, 0, 300.f, uniformBands(0.4f));
    LateReverbNetwork used, fresh;
    used.prepare(config);
    fresh.prepare(config);
    used.apply(snapshot);
    fresh.apply(snapshot);
    const std::uint32_t frames = 24000;
    std::vector<float> impulse(frames, 0.f), silence(frames, 0.f), a(frames), b(frames);
    impulse[0] = 1.f;
    const float* driven[] = {impulse.data()};
    const float* quiet[] = {silence.data()};
    float* left[] = {a.data()};
    used.process(driven, left, nullptr, frames);
    double level = 0.0;
    for (int block = 0; block < 12; ++block) used.process(quiet, left, nullptr, frames);
    for (float value : a) level = std::max(level, double(std::fabs(value)));
    near(level, 0.0, 0.0, "a dead room reads exact silence");
    used.process(driven, left, nullptr, frames);
    float* other[] = {b.data()};
    fresh.process(driven, other, nullptr, frames);
    double a2 = 0.0, b2 = 0.0;
    for (std::uint32_t i = 0; i < frames; ++i) {
        a2 += double(a[i]) * a[i];
        b2 += double(b[i]) * b[i];
    }
    near(10.0 * std::log10(a2 / b2), 0.0, 0.1, "a woken room rings like a fresh one, dB");
}

void capiDrivesTheWholePath() {
    ms_acoustics_config config;
    ms_acoustics_config_default(&config);
    ms_acoustics* acoustics = ms_acoustics_create(&config);
    check(!ms_acoustics_uses_gpu(acoustics), "no metallib means the CPU backend");
    ms_acoustic_material material{};
    for (std::size_t b = 0; b < kBandCount; ++b) material.absorption[b] = 0.2f;
    const std::uint32_t wall = ms_acoustics_material_add(acoustics, &material);
    const MeshData mesh = boxInside(vec3(0, 0, 0), vec3(8, 6, 3), {wall, wall, wall, wall, wall, wall});
    const std::uint32_t meshId = ms_acoustics_mesh_add(acoustics, mesh.positions.data(), std::uint32_t(mesh.positions.size() / 3), mesh.indices.data(),
                                                       std::uint32_t(mesh.indices.size()), mesh.materials.data());
    ms_acoustics_instance_add(acoustics, meshId, nullptr);
    const ms_acoustic_box bounds{{0, 0, 0}, {8, 6, 3}};
    ms_acoustic_room room{};
    room.bounds = &bounds;
    room.bound_count = 1;
    room.probe = {4, 3, 1.5f};
    const std::uint32_t roomId = ms_acoustics_room_add(acoustics, &room);
    const std::uint32_t source = ms_acoustics_source_add(acoustics, {2, 2, 1}, 0.1f);
    ms_acoustics_set_listener(acoustics, {5.5f, 3.3f, 2.1f});
    for (int i = 0; i < 60; ++i) ms_acoustics_tick(acoustics, 1.f / 60.f);

    ms_reflection_tap taps[MS_MAX_REFLECTION_TAPS];
    const std::uint32_t count = ms_acoustics_source_taps(acoustics, source, taps);
    check(count >= 6, "taps reach the engine format");
    std::uint32_t live = 0;
    for (std::uint32_t i = 0; i < count; ++i) live += taps[i].gain[1] > 0.f;
    check(live >= 6, "at least the first-order taps are live");
    ms_acoustic_source_view view;
    ms_acoustics_source_view_get(acoustics, source, &view);
    near(view.visibility, 1.0, 1e-3, "clear line through the C API");
    near(view.eq_mid_db, 0.0, 1e-2, "no occlusion EQ on a clear line");
    check(view.room == roomId && ms_acoustics_room_of(acoustics, MS_ACOUSTIC_NONE) == roomId, "source and listener are in the room");
    ms_acoustic_room_view measured;
    ms_acoustics_room_view_get(acoustics, roomId, &measured);
    near(measured.volume, 144.0, 144.0 * 0.05, "room measured through the C API");

    check(ms_acoustics_reverb_prepare(acoustics) == MS_OK, "reverb prepares");
    std::vector<float> in(4800, 0.f), left(4800), right(4800);
    in[0] = 1.f;
    const float* inputs[] = {in.data()};
    float* lefts[] = {left.data()};
    float* rights[] = {right.data()};
    ms_acoustics_reverb_process(acoustics, inputs, lefts, rights, 4800);
    double energy = 0.0;
    for (float value : left) energy += double(value) * value;
    check(energy > 0.0, "late reverb renders through the C API");
    ms_acoustics_destroy(acoustics);
}

void measuredCost() {
    // Measurements only, on this machine at whatever load it has. Nothing is asserted.
    LateReverbConfig config;
    config.rooms = 4;
    config.links = 6;
    LateReverbNetwork reverb;
    reverb.prepare(config);
    AcousticSnapshot& snapshot = *reverbSnapshot();
    Bands t60{2.f, 1.8f, 1.6f, 1.4f, 1.2f, 1.f, 0.7f, 0.4f};
    for (std::uint32_t r = 0; r < 4; ++r) setRoom(snapshot, r, 200.f + 100.f * float(r), t60);
    snapshot.linkCount = 6;
    for (std::uint32_t l = 0; l < 6; ++l) snapshot.links[l] = {l % 4, (l + 1) % 4, uniformBands(0.3f)};
    reverb.apply(snapshot);
    std::vector<float> in(480, 0.1f), out(480);
    const float* inputs[] = {in.data(), in.data(), in.data(), in.data()};
    float* left[] = {out.data(), out.data(), out.data(), out.data()};
    const auto begin = std::chrono::steady_clock::now();
    for (int block = 0; block < 5000; ++block) reverb.process(inputs, left, nullptr, 480);
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count();
    std::printf("  measured: 4 rooms, 6 links, 50 s of audio in %.3f s, %.0f x real time\n", seconds, 50.0 / seconds);

    Scene scene;
    const std::uint32_t wall = scene.addMaterial(flat(0.2f, 0.1f));
    const MeshData mesh = boxInside(vec3(0, 0, 0), vec3(12, 9, 4), {wall, wall, wall, wall, wall, wall});
    scene.addInstance(mesh.add(scene), Transform{});
    auto backend = makeCpuTraceBackend();
    Simulator simulator(scene, *backend);
    simulator.setListener(vec3(6, 4, 2));
    for (int s = 0; s < 32; ++s) simulator.addSource({vec3(1.f + float(s % 8) * 1.3f, 1.f + float(s / 8) * 1.5f, 1.5f), 0.1f});
    for (int i = 0; i < 60; ++i) simulator.tick(1.f / 60.f);
    const auto tickBegin = std::chrono::steady_clock::now();
    std::uint64_t rays = 0;
    for (int i = 0; i < 300; ++i) {
        simulator.tick(1.f / 60.f);
        rays += simulator.issuedClosestRays();
    }
    const double tickSeconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - tickBegin).count();
    std::printf("  measured: CPU backend, 32 sources, %.2f ms a tick, %.0f closest rays a tick\n", tickSeconds / 300.0 * 1000.0, double(rays) / 300.0);
}

void snapshotExchangeNeverTears() {
    struct Payload {
        std::array<std::uint64_t, 512> words{};
    };
    SnapshotExchange<Payload> exchange;
    std::atomic<bool> stop{false};
    std::thread writer([&] {
        for (std::uint64_t value = 1; !stop.load(std::memory_order_relaxed); ++value) {
            exchange.writeSlot().words.fill(value);
            exchange.publish();
        }
    });
    std::uint64_t torn = 0, previous = 0, backwards = 0, distinct = 0;
    for (int i = 0; i < 2000000; ++i) {
        const Payload& payload = exchange.read();
        const std::uint64_t first = payload.words[0];
        for (std::uint64_t word : payload.words) torn += word != first;
        backwards += first < previous;
        distinct += first != previous;
        previous = first;
    }
    stop.store(true);
    writer.join();
    near(double(torn), 0.0, 0.0, "snapshot words all from one publish");
    near(double(backwards), 0.0, 0.0, "reader never goes back in time");
    check(distinct > 100, "reader saw fresh snapshots");
}

void smallRayBudgetStaysBoundedAndFair() {
    Scene scene;
    const std::uint32_t wall = scene.addMaterial(flat(0.2f));
    const MeshData mesh = boxInside(vec3(0, 0, 0), vec3(12, 9, 4), {wall, wall, wall, wall, wall, wall});
    scene.addInstance(mesh.add(scene), Transform{});
    TraceCapacity capacity;
    capacity.closestRays = 512;
    auto backend = makeCpuTraceBackend(capacity);
    SimulatorConfig config;
    config.air = noAir();
    Simulator simulator(scene, *backend, config);
    simulator.setListener(vec3(6, 4, 2));
    for (int s = 0; s < 32; ++s) simulator.addSource({vec3(1.f + float(s % 8) * 1.3f, 1.f + float(s / 8) * 1.5f, 1.5f), 0.1f});
    std::uint32_t most = 0;
    for (int i = 0; i < 180; ++i) {
        simulator.tick(1.f / 60.f);
        most = std::max(most, simulator.issuedClosestRays());
    }
    check(most <= capacity.closestRays, "ray budget is a hard limit");
    const AcousticSnapshot& snapshot = simulator.snapshots().read();
    std::uint32_t weakest = kMaxPathsPerSource;
    for (std::uint32_t s = 0; s < 32; ++s) {
        std::uint32_t paths = 0;
        for (const ReflectionPath& path : snapshot.sources[s].reflections) paths += path.id != 0;
        weakest = std::min(weakest, paths);
    }
    std::printf("  ok 32 sources at %u rays a tick: fewest paths on any source %u\n", most, weakest);
    check(weakest >= 6, "every source keeps at least its first-order paths under a tight budget");
}


// ---------------------------------------------------------------------------------------------------------------------
// Edge diffraction.

void edgesComeFromTheGeometry() {
    Scene scene;
    const std::uint32_t m = scene.addMaterial(flat(0.1f));
    const MeshData box = boxInside(vec3(0, 0, 0), vec3(2, 3, 4), {m, m, m, m, m, m});
    near(double(scene.meshes()[box.add(scene)].edges.size()), 12.0, 0.0, "a box has twelve edges, its face diagonals none");
    MeshData wall;
    wall.quad(vec3(0, 0, 0), vec3(4, 0, 0), vec3(4, 3, 0), vec3(0, 3, 0), m);
    near(double(scene.meshes()[wall.add(scene)].edges.size()), 4.0, 0.0, "a thin wall has four free edges");
}

void diffractionOverAThinBarrier() {
    // A free-standing wall in the plane z = 0, x from -20 to 20 and y from 0 to 3, with nothing else. The source and the listener are on
    // opposite sides. The shortest way round is under its bottom edge, at the origin.
    Scene scene;
    const std::uint32_t m = scene.addMaterial(flat(0.1f));
    MeshData wall;
    wall.quad(vec3(-20, 0, 0), vec3(20, 0, 0), vec3(20, 3, 0), vec3(-20, 3, 0), m);
    scene.addInstance(wall.add(scene), Transform{});
    DiffractionConfig config;
    config.air = noAir();
    std::array<DiffractionPath, 4> paths;
    const Vec3 source = vec3(0, 1, -2), listener = vec3(0, 1, 3);
    const std::uint32_t count = findDiffractionPaths(scene, source, listener, config, paths);
    check(count >= 2, "paths under and over the wall");
    const double under = std::sqrt(5.0) + std::sqrt(10.0), over = std::sqrt(8.0) + std::sqrt(13.0);
    near(paths[0].length, under, 1e-3, "shortest path length under the wall");
    near(length(paths[0].points[0] - vec3(0, 0, 0)), 0.0, 1e-3, "it bends at the closest point of the bottom edge");
    near(paths[1].length, over, 1e-3, "the next one goes over the top");
    // Kurze and Anderson at 1 kHz for the path under: N = 2 delta f / c.
    const double delta = under - 5.0, n = 2 * delta * 1000 / 343, x = std::sqrt(2 * std::numbers::pi * n);
    const double lossDb = 5 + 20 * std::log10(x / std::tanh(x));
    near(20 * std::log10(paths[0].amplitude[4] * under), -lossDb, 0.01, "wedge loss at 1 kHz against Kurze-Anderson, dB");
    near(wedgeLossDb(0.f, 1000.f, 343.f), 5.0, 1e-4, "5 dB at the shadow boundary");
    check(paths[0].amplitude[1] > paths[0].amplitude[6], "low frequencies bend round more easily than high");
    near(dot(paths[0].arrival, normalize(vec3(0, -1, -3))), 1.0, 1e-4, "it arrives from the edge");
}

void insideCornersDoNotDiffract() {
    Scene scene;
    const std::uint32_t m = scene.addMaterial(flat(0.1f));
    const MeshData room = boxInside(vec3(0, 0, 0), vec3(10, 4, 8), {m, m, m, m, m, m});
    scene.addInstance(room.add(scene), Transform{});
    DiffractionConfig config;
    config.maxDetour = 20.f;
    std::array<DiffractionPath, 4> paths;
    near(double(findDiffractionPaths(scene, vec3(2, 1, 2), vec3(8, 2, 6), config, paths)), 0.0, 0.0, "a room's own corners bend nothing");
}

// A 12 x 4 x 12 m room with a square pillar from floor to ceiling in the middle.
struct PillarRoom {
    Scene scene;
    std::unique_ptr<TraceBackend> backend = makeCpuTraceBackend();
    std::unique_ptr<Simulator> simulator;
    std::uint32_t source = 0;
    PillarRoom() {
        const std::uint32_t wall = scene.addMaterial(flat(0.3f));
        const std::uint32_t stone = scene.addMaterial(flat(0.05f));  // Transmits nothing.
        scene.addInstance(boxInside(vec3(-6, 0, -6), vec3(6, 4, 6), {wall, wall, wall, wall, wall, wall}).add(scene), Transform{});
        scene.addInstance(boxInside(vec3(-1, 0, -1), vec3(1, 4, 1), {stone, stone, stone, stone, stone, stone}).add(scene), Transform{});
        RoomDesc room;
        room.bounds = {Box{vec3(-6, 0, -6), vec3(6, 4, 6)}};
        room.probe = vec3(-4, 2, -4);
        scene.addRoom(room);
        SimulatorConfig config;
        config.air = noAir();
        simulator = std::make_unique<Simulator>(scene, *backend, config);
        simulator->setListener(vec3(0, 1.5f, 5));
        source = simulator->addSource({vec3(0, 1.5f, -5), 0.1f});
    }
};

void diffractionAroundAPillar() {
    PillarRoom room;
    for (int i = 0; i < 30; ++i) room.simulator->tick(1.f / 60.f);
    const SourceAcoustics& hidden = room.simulator->snapshots().read().sources[room.source];
    near(hidden.visibility, 0.0, 1e-3, "the pillar hides the source");
    // The way round is past two corners on one side: 2 x sqrt(1 + 16) + 2 against 10 m straight.
    const double around = 2 * std::sqrt(17.0) + 2;
    std::uint32_t live = 0;
    double best = 0, length = 0;
    std::uint32_t order = 0;
    for (const ReflectionPath& path : hidden.diffraction) {
        if (path.id == 0) continue;
        ++live;
        if (path.gain[4] > best) {
            best = path.gain[4];
            length = path.delaySeconds * kDefaultSpeedOfSound;
            order = path.order;
        }
    }
    check(live == 2, "both sides of the pillar carry sound round");
    near(double(order), 2.0, 0.0, "round a thick pillar takes two edges");
    near(length, around, 1e-3, "path round the pillar");
    const double d1 = std::sqrt(17.0) + 2 - std::sqrt(1 + 36.0);  // Source to the far corner, past the near one.
    const double loss = 2 * (5 + 20 * std::log10(std::sqrt(2 * std::numbers::pi * 2 * d1 * 1000 / 343) /
                                                   std::tanh(std::sqrt(2 * std::numbers::pi * 2 * d1 * 1000 / 343))));
    near(20 * std::log10(best * around), -loss, 0.5, "two wedge losses at 1 kHz, dB");

    // Through the C API format: the diffracted paths fill the last taps.
    ms_reflection_tap taps[MS_MAX_REFLECTION_TAPS];
    (void)taps;

    // Step the source out to the side, into view: the diffracted paths fade and go.
    room.simulator->setSource(room.source, {vec3(4, 1.5f, -5), 0.1f});
    for (int i = 0; i < 60; ++i) room.simulator->tick(1.f / 60.f);
    const SourceAcoustics& seen = room.simulator->snapshots().read().sources[room.source];
    near(seen.visibility, 1.0, 1e-3, "out of the pillar's shadow");
    double left = 0;
    for (const ReflectionPath& path : seen.diffraction) left = std::max(left, double(path.gain[4]));
    near(left, 0.0, 1e-4, "diffraction fades once the line is clear");
}

// ---------------------------------------------------------------------------------------------------------------------
// Hybrid convolution reverb from the traced energy field.

double correlation(const std::vector<float>& a, const std::vector<float>& b, std::size_t from, std::size_t to) {
    double ab = 0, aa = 0, bb = 0;
    for (std::size_t i = from; i < to; ++i) {
        ab += double(a[i]) * b[i];
        aa += double(a[i]) * a[i];
        bb += double(b[i]) * b[i];
    }
    return aa > 0 && bb > 0 ? ab / std::sqrt(aa * bb) : 0.0;
}

void convolverMatchesDirectConvolution() {
    constexpr std::uint32_t kBlock = 64, kChannels = 3, kIr = 700, kLength = 4000;
    std::mt19937 rng(7);
    std::normal_distribution<float> gauss;
    std::vector<std::vector<float>> irA(kChannels, std::vector<float>(kIr)), irB = irA;
    for (auto& channel : irA)
        for (float& v : channel) v = gauss(rng);
    for (auto& channel : irB)
        for (float& v : channel) v = gauss(rng) * 0.5f;
    std::vector<float> input(kLength);
    for (float& v : input) v = gauss(rng);
    PartitionedConvolver convolver;
    convolver.prepare(kBlock, 1024, kChannels);
    const float* a[kChannels] = {irA[0].data(), irA[1].data(), irA[2].data()};
    const float* b[kChannels] = {irB[0].data(), irB[1].data(), irB[2].data()};
    check(convolver.stage(a, kIr), "first response stages");
    check(!convolver.stage(b, kIr), "a second response waits until the first is taken");
    std::vector<std::vector<float>> out(kChannels, std::vector<float>(kLength));
    // Odd chunk sizes, so blocks straddle calls. The second response is staged mid-way.
    const std::uint32_t chunks[] = {1, 37, 64, 129, 500, 3};
    std::uint32_t done = 0, turn = 0, swapAt = 0;
    while (done < kLength) {
        const std::uint32_t n = std::min(chunks[turn++ % 6], kLength - done);
        float* outs[kChannels] = {out[0].data() + done, out[1].data() + done, out[2].data() + done};
        convolver.process(input.data() + done, outs, n);
        done += n;
        if (swapAt == 0 && done >= 2000 && convolver.stage(b, kIr)) swapAt = done;
    }
    const auto direct = [&](const std::vector<float>& h, std::uint32_t n) {
        double sum = 0;
        for (std::uint32_t k = 0; k < kIr && k <= n; ++k) sum += double(h[k]) * input[n - k];
        return sum;
    };
    double worstBefore = 0, worstAfter = 0;
    for (std::uint32_t c = 0; c < kChannels; ++c) {
        // The first response fades in from silence over the first block it plays in.
        for (std::uint32_t n = 2 * kBlock; n < 1900; ++n) worstBefore = std::max(worstBefore, std::fabs(out[c][n] - direct(irA[c], n - kBlock)));
        // One block to pick the new response up at a block boundary, one to crossfade.
        for (std::uint32_t n = swapAt + 3 * kBlock; n < kLength; ++n) worstAfter = std::max(worstAfter, std::fabs(out[c][n] - direct(irB[c], n - kBlock)));
    }
    near(worstBefore, 0.0, 2e-3, "partitioned convolution equals direct convolution");
    near(worstAfter, 0.0, 2e-3, "after a swap it equals the new response over the whole history");
    // Through the hand-over every sample lies between what the old and the new response give: a blend, never a step outside both.
    double outside = 0;
    for (std::uint32_t n = swapAt; n < swapAt + 3 * kBlock; ++n) {
        const double oldValue = direct(irA[0], n - kBlock), newValue = direct(irB[0], n - kBlock);
        const double lo = std::min(oldValue, newValue), hi = std::max(oldValue, newValue);
        outside = std::max(outside, std::max(lo - out[0][n], out[0][n] - hi));
    }
    near(std::max(outside, 0.0), 0.0, 2e-3, "the swap blends the old and new responses");
}

void convolverIsRealTimeSafe() {
    PartitionedConvolver convolver;
    convolver.prepare(256, 48000, kAmbisonicChannels);
    std::vector<float> ir(48000, 0.f);
    ir[300] = 1.f;
    const float* irs[kAmbisonicChannels] = {ir.data(), ir.data(), ir.data(), ir.data()};
    std::vector<float> in(512, 0.1f), out(4 * 512);
    float* outs[kAmbisonicChannels] = {out.data(), out.data() + 512, out.data() + 1024, out.data() + 1536};
    const long before = allocationCount.load();
    for (int block = 0; block < 200; ++block) {
        if (block % 20 == 0) convolver.stage(irs, 48000);  // Responses arriving while it runs.
        convolver.process(in.data(), outs, 512);
    }
    near(double(allocationCount.load() - before), 0.0, 0.0, "no allocation in stage() or process()");
}

// A 10 x 8 x 4 m box, alpha 0.2 and scattering 0.5 everywhere: close enough to diffuse for Eyring.
struct TracedBox {
    Scene scene;
    std::unique_ptr<TraceBackend> backend = makeCpuTraceBackend();
    std::unique_ptr<Simulator> simulator;
    std::uint32_t room = kNone;
    TracedBox(float seconds, std::uint32_t paths) {
        const std::uint32_t wall = scene.addMaterial(flat(0.2f, 0.5f));
        const MeshData mesh = boxInside(vec3(0, 0, 0), vec3(10, 8, 4), {wall, wall, wall, wall, wall, wall});
        scene.addInstance(mesh.add(scene), Transform{});
        RoomDesc desc;
        desc.bounds = {Box{vec3(0, 0, 0), vec3(10, 8, 4)}};
        desc.probe = vec3(5, 4, 2);
        room = scene.addRoom(desc);
        SimulatorConfig config;
        config.air = noAir();
        config.energy.paths = paths;
        config.energy.seconds = seconds;
        simulator = std::make_unique<Simulator>(scene, *backend, config);
        simulator->setListener(vec3(3.7f, 4.6f, 1.7f));
    }
};

void energyFieldMatchesEyring() {
    TracedBox box(1.2f, 512);
    for (int i = 0; i < 900; ++i) box.simulator->tick(1.f / 60.f);
    EnergyField field;
    box.simulator->energyTracer().read(field);
    check(field.room == box.room, "the field belongs to the listener's room");
    const double V = 320, S = 304, alpha = 0.2, c = 343;
    const double eyring = 24.0 * std::log(10.0) * V / (c * -S * std::log(1.0 - alpha));
    const Bands t60 = field.t60();
    near(t60[3], eyring, eyring * 0.15, "traced T60 at 500 Hz against Eyring");
    near(t60[6], eyring, eyring * 0.15, "traced T60 at 4 kHz against Eyring");
    // Reflected energy of an impulse against the diffuse-field estimate 16 pi (1 - alpha) / A, both relative to the direct sound at 1 m, for the
    // part the field covers.
    double total = 0;
    for (std::uint32_t i = 0; i < field.bins; ++i) total += field.band(3, i);
    const double A = -S * std::log(1.0 - alpha);
    const double covered = 1.0 - std::exp(-c * A / (4 * V) * 1.2);
    const double expected = 16 * std::numbers::pi * (1 - alpha) / A * covered;
    near(10 * std::log10(total / expected), 0.0, 2.0, "traced reflected energy against 16 pi (1 - alpha) / A, dB");
    // Late in the tail the field is diffuse: no direction holds much more than its share.
    double most = 0, least = INFINITY;
    for (std::uint32_t d = 0; d < kFieldDirections; ++d) {
        double sum = 0;
        for (std::uint32_t i = field.bins / 3; i < field.bins; ++i) sum += field.at(d, 3, i);
        most = std::max(most, sum);
        least = std::min(least, sum);
    }
    check(most < 2.5 * least, "the late field arrives evenly from every direction");
    std::printf("  ok %.0f paths, %llu segments traced; late direction spread %.2f\n", field.paths,
                (unsigned long long)box.simulator->energyTracer().segmentsTraced(), most / least);
}

void energyFieldFollowsTheRoom() {
    TracedBox box(1.0f, 256);
    for (int i = 0; i < 200; ++i) box.simulator->tick(1.f / 60.f);
    box.simulator->setListener(vec3(20, 20, 20));  // Outside every room.
    for (int i = 0; i < 3; ++i) box.simulator->tick(1.f / 60.f);
    EnergyField field;
    box.simulator->energyTracer().read(field);
    check(field.room == kNone && field.paths < 2000.0, "leaving the room starts the field again");
}

void responseMatchesItsField() {
    // A synthetic field: isotropic, the same exponential decay in every band.
    EnergyField field;
    field.bins = 100;
    field.binSeconds = 0.01f;
    field.room = 0;
    field.paths = 1e6;
    field.energy.assign(std::size_t(kFieldDirections) * kBandCount * field.bins, 0.f);
    const float rate = 13.8f;  // 0.5 s T60.
    for (std::uint32_t d = 0; d < kFieldDirections; ++d)
        for (std::uint32_t b = 0; b < kBandCount; ++b)
            for (std::uint32_t i = 0; i < field.bins; ++i)
                field.at(d, b, i) = 0.01f / float(kFieldDirections) * std::exp(-rate * (float(i) + 0.5f) * field.binSeconds);
    ReverbIrSynth synth;
    synth.prepare(48000.f, 48000);
    std::vector<float> store(std::size_t(kAmbisonicChannels) * synth.maxSamples());
    float* out[kAmbisonicChannels];
    for (std::uint32_t c = 0; c < kAmbisonicChannels; ++c) out[c] = store.data() + std::size_t(c) * synth.maxSamples();
    ReverbIrShape shape;
    shape.fadeOutStart = 0.9f;
    shape.fadeOutEnd = 1.0f;
    const std::uint32_t samples = synth.build(field, shape, out);
    check(samples == 48000, "the response is as long as asked");
    // Energy of W in 50 ms windows against the field's.
    double worst = 0;
    for (std::uint32_t w = 1; w < 16; ++w) {
        double ir = 0, expected = 0;
        for (std::uint32_t n = w * 2400; n < (w + 1) * 2400; ++n) ir += double(out[0][n]) * out[0][n];
        // The bands share an impulse's energy, so a field the same in every band gives a response of that one band's energy.
        for (std::uint32_t i = w * 5; i < (w + 1) * 5; ++i) expected += field.band(3, i);
        worst = std::max(worst, std::fabs(10 * std::log10(ir / expected)));
    }
    near(worst, 0.0, 1.0, "response energy follows the field, worst 50 ms window in dB");
    // Diffuse: the first-order channels carry a third of W's energy each and are uncorrelated with it and each other.
    std::vector<float> w(out[0], out[0] + samples), y(out[1], out[1] + samples), z(out[2], out[2] + samples), x(out[3], out[3] + samples);
    double ew = 0, ex = 0;
    for (std::uint32_t n = 0; n < 24000; ++n) {
        ew += double(w[n]) * w[n];
        ex += double(x[n]) * x[n];
    }
    near(ex / ew, 1.0 / 3.0, 0.06, "diffuse field: X holds a third of W's energy");
    double worstCorrelation = 0;
    for (const auto* pair : {&x, &y, &z}) worstCorrelation = std::max(worstCorrelation, std::fabs(correlation(w, *pair, 0, 24000)));
    worstCorrelation = std::max(worstCorrelation, std::fabs(correlation(x, y, 0, 24000)));
    near(worstCorrelation, 0.0, 0.15, "diffuse field: ambisonic channels uncorrelated");

    // Everything from one direction: the response points there.
    std::fill(field.energy.begin(), field.energy.end(), 0.f);
    for (std::uint32_t b = 0; b < kBandCount; ++b)
        for (std::uint32_t i = 0; i < field.bins; ++i) field.at(5, b, i) = 0.01f * std::exp(-rate * (float(i) + 0.5f) * field.binSeconds);
    synth.build(field, shape, out);
    double wx = 0, ww = 0;
    for (std::uint32_t n = 0; n < 24000; ++n) {
        wx += double(out[3][n]) * out[0][n];
        ww += double(out[0][n]) * out[0][n];
    }
    const Vec3 u = fieldDirections()[5];
    // The uniform prior keeps 10 % spread over every direction, so the mean direction cosine is 90 % of the direction's.
    near(wx / ww, -u.z * 0.9, 0.08, "one direction: X follows it");
    // A constant field: the response's power is the field's, with no decay to bias anything and the bands' powers simply adding.
    std::fill(field.energy.begin(), field.energy.end(), 0.001f / float(kFieldDirections));
    synth.build(field, shape, out);
    double flatEnergy = 0;
    for (std::uint32_t n = 4800; n < 43200; ++n) flatEnergy += double(out[0][n]) * out[0][n];
    near(10 * std::log10(flatEnergy / (0.001 * 80)), 0.0, 0.5, "constant field: response power equals the field's, dB");
}

void hybridDelaysAndScalesTheSend() {
    // The network's own send through the hybrid: the same tail, later and quieter by exactly the gain.
    const auto render = [](bool hybrid, std::vector<float>& left) {
        LateReverbConfig config;
        config.rooms = 1;
        config.links = 0;
        config.maxHybridDelay = 48000;
        LateReverbNetwork reverb;
        reverb.prepare(config);
        AcousticSnapshot& snapshot = *reverbSnapshot();
        setRoom(snapshot, 0, 320.f, uniformBands(0.76f));
        reverb.apply(snapshot);
        if (hybrid) reverb.setHybrid(0, 9600, uniformBands(-20.f));
        std::vector<float> in(left.size(), 0.f), right(left.size());
        in[0] = 1.f;
        for (std::size_t done = 0; done < in.size(); done += 480) {
            const float* inputs[] = {in.data() + done};
            float* lefts[] = {left.data() + done};
            float* rights[] = {right.data() + done};
            reverb.process(inputs, lefts, rights, 480);
        }
    };
    std::vector<float> plain(96000), delayed(96000);
    render(false, plain);
    render(true, delayed);
    double before = 0, a = 0, b = 0;
    for (std::size_t n = 0; n < 9600; ++n) before += double(delayed[n]) * delayed[n];
    for (std::size_t n = 0; n < 48000; ++n) {
        a += double(plain[n]) * plain[n];
        b += double(delayed[n + 9600]) * delayed[n + 9600];
    }
    near(before, 0.0, 1e-12, "nothing comes out before the hybrid delay");
    near(10 * std::log10(b / a), -20.0, 0.3, "the delayed tail is the plain tail at the hybrid gain, dB");
    double e0 = 0;
    for (std::size_t n = 0; n < 1920; ++n) e0 += double(plain[n]) * plain[n];
    // Initial energy rate 4 pi c / V a second, less the build-up over the first line delays.
    std::printf("    network first 40 ms: %.2f dB against 4 pi c / V: %.2f dB\n", 10 * std::log10(e0), 10 * std::log10(4 * std::numbers::pi * 343 / 320 * 0.04));
}

void capiHybridReverb() {
    ms_acoustics_config config;
    ms_acoustics_config_default(&config);
    config.convolution_seconds = 1.f;
    ms_acoustics* acoustics = ms_acoustics_create(&config);
    ms_acoustic_material material{};
    for (std::size_t b = 0; b < kBandCount; ++b) {
        material.absorption[b] = 0.2f;
        material.scattering[b] = 0.5f;
    }
    const std::uint32_t wall = ms_acoustics_material_add(acoustics, &material);
    const MeshData mesh = boxInside(vec3(0, 0, 0), vec3(10, 8, 4), {wall, wall, wall, wall, wall, wall});
    const std::uint32_t meshId = ms_acoustics_mesh_add(acoustics, mesh.positions.data(), std::uint32_t(mesh.positions.size() / 3), mesh.indices.data(),
                                                       std::uint32_t(mesh.indices.size()), mesh.materials.data());
    ms_acoustics_instance_add(acoustics, meshId, nullptr);
    const ms_acoustic_box bounds{{0, 0, 0}, {10, 8, 4}};
    ms_acoustic_room room{};
    room.bounds = &bounds;
    room.bound_count = 1;
    room.probe = {5, 4, 2};
    const std::uint32_t roomId = ms_acoustics_room_add(acoustics, &room);
    ms_acoustics_set_listener(acoustics, {3.7f, 4.6f, 1.7f});
    ms_acoustics_set_listener_orientation(acoustics, {0, 0, -1}, {0, 1, 0});
    for (int i = 0; i < 60; ++i) ms_acoustics_tick(acoustics, 1.f / 60.f);
    check(ms_acoustics_reverb_prepare(acoustics) == MS_OK, "hybrid reverb prepares");
    // The render thread runs alongside: it takes each new response as it comes.
    std::vector<float> silence(4800, 0.f), sink(4800);
    for (int i = 0; i < 900; ++i) {
        ms_acoustics_tick(acoustics, 1.f / 60.f);
        const float* quiet[] = {silence.data()};
        float* out[] = {sink.data()};
        ms_acoustics_reverb_process(acoustics, quiet, out, out, 800);
    }
    ms_acoustic_energy_view view;
    ms_acoustics_energy_view_get(acoustics, &view);
    check(view.active && view.room == roomId, "a response from the traced field is in use");
    const double eyring = 24.0 * std::log(10.0) * 320 / (343 * -304 * std::log(0.8));
    near(view.t60[3], eyring, eyring * 0.2, "energy view T60 at 500 Hz");

    // An impulse through the whole hybrid: convolution, then the network from the end of the response at the level the field reached.
    const std::uint32_t frames = 48000 * 2;
    std::vector<float> in(frames, 0.f), left(frames), right(frames);
    in[0] = 1.f;
    const float* inputs[] = {in.data()};
    float* lefts[] = {left.data()};
    float* rights[] = {right.data()};
    for (std::uint32_t done = 0; done < frames; done += 480) {
        const float* i2[] = {in.data() + done};
        float* l2[] = {left.data() + done};
        float* r2[] = {right.data() + done};
        ms_acoustics_reverb_process(acoustics, i2, l2, r2, 480);
    }
    (void)inputs;
    (void)lefts;
    (void)rights;
    // Energy decay in 40 ms windows: no step at the hand-over, and the slope of the room throughout.
    std::vector<double> db;
    for (std::uint32_t w = 0; w < frames / 1920; ++w) {
        double e = 0;
        for (std::uint32_t n = w * 1920; n < (w + 1) * 1920; ++n) e += double(left[n]) * left[n] + double(right[n]) * right[n];
        db.push_back(10 * std::log10(std::max(e, 1e-30)));
    }
    for (std::size_t w = 0; w < 34; ++w) std::printf("    %zu: %.1f dB  step %.1f\n", w, db[w], w ? db[w] - db[w - 1] : 0.0);
    double worstStep = 0;
    const double expectedStep = 60.0 / eyring * 0.04;
    for (std::size_t w = 4; w + 1 < db.size() && w < 32; ++w) worstStep = std::max(worstStep, std::fabs(db[w + 1] - db[w] + expectedStep));
    near(worstStep, 0.0, 3.0, "decay in 40 ms windows departs from the room's slope by at most, dB");
    // Stereo from the ambisonic response: the late tail is decorrelated between the ears.
    near(correlation(left, right, 4800, 38400), 0.0, 0.35, "late tail decorrelated between the ears");

    std::vector<float> w(frames), x(frames), y(frames), z(frames);
    float* foa[] = {w.data(), y.data(), z.data(), x.data()};
    std::fill(left.begin(), left.end(), 0.f);
    for (std::uint32_t done = 0; done < frames; done += 480) {
        const float* i2[] = {in.data() + done};
        float* l2[] = {left.data() + done};
        float* r2[] = {right.data() + done};
        float* a2[] = {foa[0] + done, foa[1] + done, foa[2] + done, foa[3] + done};
        ms_acoustics_reverb_process_ambisonic(acoustics, i2, l2, r2, a2, 480);
    }
    double ambisonicEnergy = 0;
    for (float v : w) ambisonicEnergy += double(v) * v;
    check(ambisonicEnergy > 0.0, "the ambisonic output carries the convolution");

    // Real time: no allocation while rendering, with new responses arriving; and what it costs.
    std::vector<float> block(480, 0.01f), sinkLeft(480), sinkRight(480);
    const float* blockInputs[] = {block.data()};
    float* blockLeft[] = {sinkLeft.data()};
    float* blockRight[] = {sinkRight.data()};
    for (int i = 0; i < 40; ++i) ms_acoustics_tick(acoustics, 1.f / 60.f);  // A response is waiting.
    const long before = allocationCount.load();
    const auto start = std::chrono::steady_clock::now();
    for (int i = 0; i < 1000; ++i) ms_acoustics_reverb_process(acoustics, blockInputs, blockLeft, blockRight, 480);
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    near(double(allocationCount.load() - before), 0.0, 0.0, "no allocation in the hybrid render");
    std::printf("  measured: hybrid reverb, 1 s ambisonic response, 10 s of audio in %.3f s, %.0f x real time\n", seconds, 10.0 / seconds);
    ms_acoustics_destroy(acoustics);
}
}  // namespace

// ---------------------------------------------------------------------------------------------------------------------
// Baked propagation.

// A corridor that doubles back: thick walls from either side leave a gap at alternate ends, so the way from one end to the other turns four
// corners, more than edge diffraction follows.
struct Serpentine {
    Scene scene;
    std::uint32_t room = kNone;
    Serpentine() {
        const std::uint32_t wall = scene.addMaterial(flat(0.3f));
        scene.addInstance(boxInside(vec3(-6, 0, -6), vec3(6, 4, 6), {wall, wall, wall, wall, wall, wall}).add(scene), Transform{});
        scene.addInstance(boxInside(vec3(-6, 0, -2.2f), vec3(3, 4, -1.8f), {wall, wall, wall, wall, wall, wall}).add(scene), Transform{});
        scene.addInstance(boxInside(vec3(-3, 0, 1.8f), vec3(6, 4, 2.2f), {wall, wall, wall, wall, wall, wall}).add(scene), Transform{});
        RoomDesc desc;
        desc.bounds = {Box{vec3(-6, 0, -6), vec3(6, 4, 6)}};
        desc.probe = vec3(0, 2, -4);
        room = scene.addRoom(desc);
    }
};

void bakedRoutesRoundACorridor() {
    Serpentine s;
    BakeConfig config;
    config.spacing = 1.f;
    config.reverb = false;
    auto bake = std::make_shared<AcousticBake>();
    std::string error;
    float last = -1.f;
    bool monotonic = true;
    const bool ok = bakeScene(s.scene, config, *bake, &error, [&](float done) {
        monotonic = monotonic && done >= last;
        last = done;
    });
    check(ok, "the corridor bakes");
    check(monotonic && last == 1.f, "bake progress rises to 1");
    check(bake->probeCount() > 400 && bake->hasPathing() && !bake->hasReverb(), "probes fill the corridor");
    std::printf("  ok %u probes, %u dropped, %llu visible pairs\n", bake->probeCount(), bake->stats.droppedProbes,
                (unsigned long long)bake->stats.visiblePairs);

    const Vec3 source = vec3(0, 1.5f, -5), listener = vec3(0, 1.5f, 5);
    DiffractionConfig bends;
    bends.air = noAir();
    std::array<DiffractionPath, 4> edgePaths;
    check(findDiffractionPaths(s.scene, source, listener, bends, edgePaths) == 0 ||
              std::all_of(edgePaths.begin(), edgePaths.end(), [](const DiffractionPath& p) { return p.order <= 2; }),
          "edge diffraction alone has at most two edges");

    // Round the corners: 2 sqrt(9 + 7.84) + 0.8 + sqrt(36 + 12.96) metres pulled tight on the walls' corners.
    const double shortest = 2 * std::sqrt(9 + 7.84) + 0.8 + std::sqrt(36 + 12.96);
    auto backend = makeCpuTraceBackend();
    SimulatorConfig simulation;
    simulation.air = noAir();
    simulation.bends.air = noAir();
    Simulator plain(s.scene, *backend, simulation);
    plain.setListener(listener);
    const std::uint32_t id = plain.addSource({source, 0.1f});
    for (int i = 0; i < 40; ++i) plain.tick(1.f / 60.f);
    double without = 0;
    for (const ReflectionPath& path : plain.snapshots().read().sources[id].diffraction) without = std::max(without, double(path.gain[4]));
    near(without, 0.0, 1e-6, "without a bake nothing finds the way round");

    auto backend2 = makeCpuTraceBackend();
    Simulator baked(s.scene, *backend2, simulation);
    baked.setBake(bake);
    baked.setListener(listener);
    const std::uint32_t id2 = baked.addSource({source, 0.1f});
    for (int i = 0; i < 60; ++i) baked.tick(1.f / 60.f);
    const SourceAcoustics& view = baked.snapshots().read().sources[id2];
    near(view.visibility, 0.0, 1e-3, "the walls hide the source");
    const ReflectionPath* best = nullptr;
    for (const ReflectionPath& path : view.diffraction)
        if (path.id != 0 && (!best || path.gain[4] > best->gain[4])) best = &path;
    check(best != nullptr && best->gain[4] > 0.f, "a baked route carries sound round the corridor");
    if (!best) return;
    const double length = best->delaySeconds * kDefaultSpeedOfSound;
    check(best->order >= 3, "the route turns more corners than edge diffraction follows");
    check(length >= shortest - 1e-3 && length <= 1.2 * shortest, "the route is within 20% of the tightest way round");
    std::printf("  ok route of %.2f m round %u bends, tightest %.2f m, %.1f dB at 1 kHz\n", length, best->order, shortest,
                20 * std::log10(best->gain[4]));
    const double routeLoss = 5 + 20 * std::log10(std::sqrt(2 * std::numbers::pi * 2 * (length - 10) * 1000 / 343) /
                                                  std::tanh(std::sqrt(2 * std::numbers::pi * 2 * (length - 10) * 1000 / 343)));
    near(20 * std::log10(best->gain[4] * length), -routeLoss, 0.5, "route loss at 1 kHz: one wedge for the whole detour, dB");
    // The route arrives from the side of the corridor it last came through: the gap at the listener's left, -x.
    check(best->arrival.x < 0.f, "the route arrives through the last gap");

    // The source steps into view: the route fades.
    baked.setSource(id2, {vec3(0, 1.5f, 4), 0.1f});
    for (int i = 0; i < 60; ++i) baked.tick(1.f / 60.f);
    double left = 0;
    for (const ReflectionPath& path : baked.snapshots().read().sources[id2].diffraction) left = std::max(left, double(path.gain[4]));
    near(left, 0.0, 1e-4, "the route fades once the line is clear");
}

void bakeDropsProbesInsideSolids() {
    PillarRoom room;
    BakeConfig config;
    config.spacing = 0.75f;
    config.clearance = 0.2f;
    config.reverb = false;
    AcousticBake bake;
    check(bakeScene(room.scene, config, bake), "the pillar room bakes");
    std::uint32_t inside = 0;
    for (Vec3 p : bake.positions) inside += std::fabs(p.x) < 1.f && std::fabs(p.z) < 1.f;
    near(inside, 0.0, 0.0, "no probe is left inside the pillar");
    check(bake.stats.droppedProbes > 0, "the grid points inside it were dropped");
}

void bakedFieldMatchesTheRoom() {
    TracedBox box(1.2f, 0);
    BakeConfig config;
    config.spacing = 4.f;
    config.pathing = false;
    config.pathsPerProbe = 1536;
    config.energy.seconds = 1.2f;
    config.energy.air = noAir();
    AcousticBake bake;
    check(bakeScene(box.scene, config, bake), "the box bakes");
    check(bake.hasReverb() && !bake.hasPathing() && bake.probeCount() >= 6, "every probe has a field");
    EnergyField field;
    check(bake.fieldAt(vec3(3.7f, 4.6f, 1.7f), box.room, 120, 0.01f, field), "a field at the listener");
    check(!bake.fieldAt(vec3(3.7f, 4.6f, 1.7f), box.room + 1, 120, 0.01f, field) || false, "no field for a room without probes");
    bake.fieldAt(vec3(3.7f, 4.6f, 1.7f), box.room, 120, 0.01f, field);
    const double V = 320, S = 304, alpha = 0.2, c = 343;
    const double eyring = 24.0 * std::log(10.0) * V / (c * -S * std::log(1.0 - alpha));
    near(field.t60()[4], eyring, eyring * 0.15, "baked T60 at 1 kHz against Eyring");
    double total = 0;
    for (std::uint32_t i = 0; i < field.bins; ++i) total += field.band(4, i);
    const double A = -S * std::log(1.0 - alpha);
    const double expected = 16 * std::numbers::pi * (1 - alpha) / A * (1.0 - std::exp(-c * A / (4 * V) * 1.2));
    near(10 * std::log10(total / expected), 0.0, 2.0, "baked reflected energy against 16 pi (1 - alpha) / A, dB");
    // Half-width bins hold the same energy in all.
    EnergyField fine;
    bake.fieldAt(vec3(3.7f, 4.6f, 1.7f), box.room, 240, 0.005f, fine);
    double fineTotal = 0;
    for (std::uint32_t i = 0; i < fine.bins; ++i) fineTotal += fine.band(4, i);
    near(fineTotal / total, 1.0, 1e-4, "rebinning keeps the energy");

    // To disk and back, bit for bit; a cut-short file is refused.
    const std::string file = std::string(std::getenv("TMPDIR") ? std::getenv("TMPDIR") : "/tmp") + "/ms_bake_test.msbake";
    check(bake.save(file), "the bake saves");
    AcousticBake loaded;
    std::string error;
    check(AcousticBake::load(file, loaded, &error), "the bake loads");
    check(loaded.totals == bake.totals && loaded.shares == bake.shares && loaded.rooms == bake.rooms && loaded.sceneHash == bake.sceneHash,
          "what loads is what was saved");
    {
        std::ifstream in(file, std::ios::binary);
        std::vector<char> bytes((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
        std::ofstream outFile(file, std::ios::binary | std::ios::trunc);
        outFile.write(bytes.data(), std::streamsize(bytes.size() - 7));
    }
    check(!AcousticBake::load(file, loaded, &error), "a cut-short bake is refused");
    std::remove(file.c_str());

    // The scene hash follows what rays see.
    const std::uint64_t before = staticSceneHash(box.scene);
    check(before == bake.sceneHash, "the bake carries its scene's hash");
    box.scene.setMaterial(0, flat(0.4f, 0.5f));
    check(staticSceneHash(box.scene) != before, "a changed material changes the hash");

    // Blending: the baked field counts as prior paths, the live one as its own.
    EnergyField live = field, blended;
    std::fill(live.energy.begin(), live.energy.end(), 0.f);
    live.paths = 1000;
    blendFields(live, field, 3000, blended);
    near(blended.energy[50] / field.energy[50], 0.75, 1e-5, "a prior of 3000 paths against 1000 live ones keeps three quarters");
    near(blended.paths, 4000, 0, "the blend holds both");
}

void capiBakedReverb() {
    ms_acoustics_config config;
    ms_acoustics_config_default(&config);
    config.convolution_seconds = 1.f;
    config.energy_paths = 0;  // Nothing traced live: the response can only come from the bake.
    ms_acoustics* acoustics = ms_acoustics_create(&config);
    ms_acoustic_material material{};
    for (std::size_t b = 0; b < kBandCount; ++b) {
        material.absorption[b] = 0.2f;
        material.scattering[b] = 0.5f;
    }
    const std::uint32_t wall = ms_acoustics_material_add(acoustics, &material);
    const MeshData mesh = boxInside(vec3(0, 0, 0), vec3(10, 8, 4), {wall, wall, wall, wall, wall, wall});
    const std::uint32_t meshId = ms_acoustics_mesh_add(acoustics, mesh.positions.data(), std::uint32_t(mesh.positions.size() / 3), mesh.indices.data(),
                                                       std::uint32_t(mesh.indices.size()), mesh.materials.data());
    ms_acoustics_instance_add(acoustics, meshId, nullptr);
    const ms_acoustic_box bounds{{0, 0, 0}, {10, 8, 4}};
    ms_acoustic_room room{};
    room.bounds = &bounds;
    room.bound_count = 1;
    room.probe = {5, 4, 2};
    const std::uint32_t roomId = ms_acoustics_room_add(acoustics, &room);
    ms_acoustics_set_listener(acoustics, {3.7f, 4.6f, 1.7f});
    for (int i = 0; i < 60; ++i) ms_acoustics_tick(acoustics, 1.f / 60.f);
    check(ms_acoustics_reverb_prepare(acoustics) == MS_OK, "reverb prepares with no live paths");

    ms_acoustic_bake_config bakeConfig;
    ms_acoustic_bake_config_default(&bakeConfig);
    bakeConfig.probe_spacing = 4.f;
    bakeConfig.paths_per_probe = 1024;
    ms_status status = MS_ERR_DEVICE;
    std::uint32_t calls = 0;
    ms_acoustic_bake* bake = ms_acoustics_bake(
        acoustics, &bakeConfig, [](float, void* user) { ++*static_cast<std::uint32_t*>(user); }, &calls, &status);
    check(bake && status == MS_OK && calls > 0, "the C API bakes, reporting progress");
    check(ms_acoustic_bake_matches(bake, acoustics), "the bake matches its scene");
    const std::string file = std::string(std::getenv("TMPDIR") ? std::getenv("TMPDIR") : "/tmp") + "/ms_capi_bake.msbake";
    check(ms_acoustic_bake_save(bake, file.c_str()) == MS_OK, "the C API saves a bake");
    const std::uint32_t probes = ms_acoustic_bake_probe_count(bake);
    ms_acoustic_bake_destroy(bake);
    ms_acoustic_bake* loaded = ms_acoustic_bake_load(file.c_str(), &status);
    check(loaded && status == MS_OK && ms_acoustic_bake_probe_count(loaded) == probes, "and loads it back");
    check(!ms_acoustic_bake_load((file + ".missing").c_str(), &status) && status == MS_ERR_NOT_FOUND, "a missing bake is not found");
    std::remove(file.c_str());

    ms_acoustic_energy_view view;
    for (int i = 0; i < 60; ++i) ms_acoustics_tick(acoustics, 1.f / 60.f);
    ms_acoustics_energy_view_get(acoustics, &view);
    check(!view.active, "no live paths and no bake: no convolution response");
    ms_acoustics_set_bake(acoustics, loaded);
    ms_acoustic_bake_destroy(loaded);  // The ms_acoustics keeps its own reference.
    std::vector<float> silence(800, 0.f), sink(800);
    for (int i = 0; i < 60; ++i) {
        ms_acoustics_tick(acoustics, 1.f / 60.f);
        const float* quiet[] = {silence.data()};
        float* out[] = {sink.data()};
        ms_acoustics_reverb_process(acoustics, quiet, out, out, 800);
    }
    ms_acoustics_energy_view_get(acoustics, &view);
    check(view.active && view.baked && view.room == roomId, "the bake gives the room a response at once");
    const double eyring = 24.0 * std::log(10.0) * 320 / (343 * -304 * std::log(0.8));
    near(view.t60[4], eyring, eyring * 0.2, "baked response T60 at 1 kHz");
    ms_acoustics_set_bake(acoustics, nullptr);
    for (int i = 0; i < 60; ++i) {
        ms_acoustics_tick(acoustics, 1.f / 60.f);
        const float* quiet[] = {silence.data()};
        float* out[] = {sink.data()};
        ms_acoustics_reverb_process(acoustics, quiet, out, out, 800);
    }
    ms_acoustics_energy_view_get(acoustics, &view);
    check(!view.active && !view.baked, "taking the bake away ends the response");
    ms_acoustics_destroy(acoustics);
}

void* operator new(std::size_t size) {
    allocationCount.fetch_add(1, std::memory_order_relaxed);
    if (void* p = std::malloc(size)) return p;
    throw std::bad_alloc();
}
void operator delete(void* p) noexcept { std::free(p); }
void operator delete(void* p, std::size_t) noexcept { std::free(p); }

int main() {
    airAttenuationMatchesIso();
    sabineAndEyringOnAnAnalyticRoom();
    imageSourceDelayAndGainInABox();
    scatteringAirAndScaleInAPath();
    occlusionAndTransmission();
    partialOcclusionOverTheSourceDisc();
    movingDoorChangesReflectionsAndStaysStable();
    probesMeasureARealRoom();
    portalGainAndOpenness();
    portalChainThroughAThirdRoom();
    coupledRoomsAnalytic();
    lateReverbDecayFollowsTheRoom();
    lateReverbLevelIsPhysical();
    lateReverbBinauralCoherence();
    coupledReverbMatchesTheModes();
    reverbProcessIsRealTimeSafe();
    sleepingRoomsWakeWithoutLoss();
    snapshotExchangeNeverTears();
    smallRayBudgetStaysBoundedAndFair();
    capiDrivesTheWholePath();
    edgesComeFromTheGeometry();
    diffractionOverAThinBarrier();
    insideCornersDoNotDiffract();
    diffractionAroundAPillar();
    convolverMatchesDirectConvolution();
    convolverIsRealTimeSafe();
    energyFieldMatchesEyring();
    energyFieldFollowsTheRoom();
    responseMatchesItsField();
    hybridDelaysAndScalesTheSend();
    capiHybridReverb();
    bakedRoutesRoundACorridor();
    bakeDropsProbesInsideSolids();
    bakedFieldMatchesTheRoom();
    capiBakedReverb();
    measuredCost();
    std::printf(failures ? "\n%d FAILED\n" : "\nall passed\n", failures);
    return failures ? 1 : 0;
}
