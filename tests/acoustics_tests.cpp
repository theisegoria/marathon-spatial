// CPU-only acoustics tests against analytic rooms. Everything renders to memory; nothing here touches an audio device or the GPU.
#include "acoustics/AcousticSimulator.hpp"
#include "acoustics/LateReverb.hpp"
#include "marathon_spatial/marathon_spatial_acoustics.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdio>
#include <cstdlib>
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

}  // namespace

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
    coupledReverbMatchesTheModes();
    reverbProcessIsRealTimeSafe();
    sleepingRoomsWakeWithoutLoss();
    snapshotExchangeNeverTears();
    smallRayBudgetStaysBoundedAndFair();
    capiDrivesTheWholePath();
    measuredCost();
    std::printf(failures ? "\n%d FAILED\n" : "\nall passed\n", failures);
    return failures ? 1 : 0;
}
