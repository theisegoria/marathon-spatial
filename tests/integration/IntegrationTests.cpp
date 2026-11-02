// SPDX-License-Identifier: Apache-2.0
// The library's parts together, through the C API a host uses. Everything renders into memory: this target links neither the
// output library nor any device code, so nothing here can reach an audio device, and no Metal device is created.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

#include "RtProbe.h"
#include "marathon_spatial/marathon_spatial.h"
#include "marathon_spatial/marathon_spatial_acoustics.h"
#include "marathon_spatial/marathon_spatial_physics.h"

namespace {

int failures = 0;

void check(bool ok, const char* what) {
    std::printf("%s %s\n", ok ? "  ok" : "FAIL", what);
    if (!ok) ++failures;
}
void near(double actual, double expected, double tolerance, const char* what) {
    const bool ok = std::fabs(actual - expected) <= tolerance;
    std::printf("%s %-62s %.6g expected %.6g (+- %.3g)\n", ok ? "  ok" : "FAIL", what, actual, expected, tolerance);
    if (!ok) ++failures;
}

constexpr float kRate = 48000.0f;
constexpr std::uint32_t kBlock = 128;

// A wide floor at y = 0: the one surface, so the one path is known exactly.
struct Floor {
    ms_acoustics* acoustics;
    std::uint32_t source;
    Floor(ms_vec3 sourcePosition, ms_vec3 listener) {
        ms_acoustics_config config;
        ms_acoustics_config_default(&config);
        acoustics = ms_acoustics_create(&config);
        ms_acoustic_material material{};
        for (float& a : material.absorption) a = 0.1f;
        const std::uint32_t floorMaterial = ms_acoustics_material_add(acoustics, &material);
        const float positions[] = {-500, 0, -500, -500, 0, 500, 500, 0, 500, 500, 0, -500};
        const std::uint32_t indices[] = {0, 1, 2, 0, 2, 3};
        const std::uint32_t materials[] = {floorMaterial, floorMaterial};
        const std::uint32_t mesh = ms_acoustics_mesh_add(acoustics, positions, 4, indices, 6, materials);
        ms_acoustics_instance_add(acoustics, mesh, nullptr);
        source = ms_acoustics_source_add(acoustics, sourcePosition, 0.1f);
        ms_acoustics_set_listener(acoustics, listener);
        for (int i = 0; i < 120; ++i) ms_acoustics_tick(acoustics, 1.0f / 60.0f);
    }
    ~Floor() { ms_acoustics_destroy(acoustics); }
    // The tap that comes up from the floor, or a tap with a negative delay if there is none.
    ms_reflection_tap floorTap() const {
        ms_reflection_tap taps[MS_MAX_REFLECTION_TAPS];
        const std::uint32_t count = ms_acoustics_source_taps(acoustics, source, taps);
        for (std::uint32_t i = 0; i < count; ++i)
            if (taps[i].gain[1] > 0.0f && taps[i].direction.y < -0.1f) return taps[i];
        ms_reflection_tap none{};
        none.delay_seconds = -1.0f;
        return none;
    }
};

double onsetSeconds(const std::vector<float>& x) {
    float peak = 0.0f;
    for (float v : x) peak = std::max(peak, std::fabs(v));
    for (std::size_t i = 0; i < x.size(); ++i)
        if (peak > 0.0f && std::fabs(x[i]) >= 0.1f * peak) return static_cast<double>(i) / kRate;
    return -1.0;
}

std::vector<float> renderImpulse(ms_vec3 position, bool propagation, const ms_reflection_tap* tap) {
    ms_engine_desc d;
    ms_engine_desc_default(&d);
    d.smoothing_seconds = 0.0f;
    ms_engine* e = nullptr;
    ms_engine_create(&d, &e);
    ms_hrtf* sphere = nullptr;
    ms_hrtf_create_sphere(0.0875f, kRate, &sphere);
    ms_engine_set_hrtf(e, sphere);
    ms_hrtf_release(sphere);
    std::vector<float> impulse(64, 0.0f);
    impulse[0] = 1.0f;
    ms_buffer b;
    ms_source s;
    ms_buffer_create(e, impulse.data(), impulse.size(), 1, kRate, &b);
    ms_source_create(e, &s);
    ms_source_params p;
    ms_source_params_default(&p);
    p.position = position;
    p.doppler = 0.0f;
    p.air_absorption = 0.0f;
    p.near_field = 0.0f;
    p.propagation_delay = propagation;
    if (tap) {  // Only the reflection: the distance curve silences the direct sound, which leaves the source's own timing alone.
        p.distance_model = MS_DISTANCE_LINEAR;
        p.ref_distance = 1.0f;
        p.max_distance = 2.0f;
        ms_reflection_tap loud = *tap;
        for (float& g : loud.gain) g = 1.0f;
        ms_source_set_reflections(e, s, &loud, 1);
    } else {
        p.distance_model = MS_DISTANCE_NONE;
    }
    ms_source_set_params(e, s, &p);
    ms_source_set_buffer(e, s, b, 0, 0);
    ms_source_play(e, s, 0);
    ms_engine_commit(e);
    std::vector<float> left(19200), right(19200);
    ms_engine_render(e, left.data(), right.data(), 19200);
    ms_engine_destroy(e);
    for (std::size_t i = 0; i < left.size(); ++i) left[i] += right[i];
    return left;
}

// Review finding 5: a tap delay is path minus direct over c, in the simulator's output and at the ears.
void tapDelayIsTheGap() {
    std::printf("== tap delay is the gap to the direct sound\n");
    const ms_vec3 source{0, 10, -34}, listener{0, 10, 0}, relative{0, 0, -34};
    const double direct = 34.0, reflected = std::sqrt(34.0 * 34.0 + 20.0 * 20.0), expected = (reflected - direct) / 343.0;
    const Floor floor(source, listener);
    const ms_reflection_tap tap = floor.floorTap();
    near(tap.delay_seconds, expected, 2e-4, "the simulator's tap at 34 m");
    for (const bool propagation : {false, true}) {
        const double directAt = onsetSeconds(renderImpulse(relative, propagation, nullptr));
        const double tapAt = onsetSeconds(renderImpulse(relative, propagation, &tap));
        near(tapAt - directAt, expected, 1e-3, propagation ? "gap at the ears, propagation delay on" : "gap at the ears, propagation delay off");
    }

    // Filtered on the same value: over 0.17 s of path, but close behind the direct sound.
    const ms_vec3 far{0, 10, -70};
    const Floor distant(far, listener);
    const ms_reflection_tap late = distant.floorTap();
    const double gap = (std::sqrt(70.0 * 70.0 + 20.0 * 20.0) - 70.0) / 343.0;
    near(late.delay_seconds, gap, 2e-4, "a path over 0.17 s long survives on its gap");

    ms_acoustic_portal_path paths[MS_ACOUSTIC_MAX_PORTAL_PATHS];
    ms_acoustics_source_portal_paths(floor.acoustics, floor.source, paths);
    check(paths[0].portals == 0, "no portals in this scene, so no portal path");
}

// Portal paths use the same definition, measured from the straight line.
void portalDelayIsTheGap() {
    std::printf("== portal path delay is the gap to the straight line\n");
    ms_acoustics_config config;
    ms_acoustics_config_default(&config);
    ms_acoustics* a = ms_acoustics_create(&config);
    const ms_acoustic_box left{{-20, 0, -5}, {0, 5, 5}}, right{{0, 0, -5}, {20, 5, 5}};
    ms_acoustic_room room{};
    room.bounds = &left;
    room.bound_count = 1;
    room.probe = {-10, 2, 0};
    const std::uint32_t roomA = ms_acoustics_room_add(a, &room);
    room.bounds = &right;
    room.probe = {10, 2, 0};
    const std::uint32_t roomB = ms_acoustics_room_add(a, &room);
    ms_acoustic_portal portal{};
    portal.room_a = roomA;
    portal.room_b = roomB;
    portal.center = {0, 2, 3};
    portal.normal = {1, 0, 0};
    portal.half_u = {0, 0, 1};
    portal.half_v = {0, 1, 0};
    ms_acoustics_portal_add(a, &portal);
    const std::uint32_t source = ms_acoustics_source_add(a, {-10, 2, 0}, 0.1f);
    ms_acoustics_set_listener(a, {10, 2, 0});
    for (int i = 0; i < 120; ++i) ms_acoustics_tick(a, 1.0f / 60.0f);
    ms_acoustic_portal_path paths[MS_ACOUSTIC_MAX_PORTAL_PATHS];
    ms_acoustics_source_portal_paths(a, source, paths);
    check(paths[0].portals == 1, "one portal path");
    // The route bends through the aperture, whose nearest edge is 2 m off the line.
    const double straight = 20.0 / 343.0;
    check(paths[0].delay_seconds >= 0.0f && paths[0].delay_seconds < 0.5 * straight, "the delay is the extra route, not the whole route");
    ms_acoustics_destroy(a);
}

// Review finding 18: the physics runs from the engine's clock whatever the sources do, and everything on the render path is clean.
void physicsOnTheBlockClock() {
    std::printf("== physics on the engine's block clock\n");
    ms_engine_desc d;
    ms_engine_desc_default(&d);
    ms_engine* e = nullptr;
    ms_engine_create(&d, &e);
    ms_hrtf* sphere = nullptr;
    ms_hrtf_create_sphere(0.0875f, kRate, &sphere);
    ms_engine_set_hrtf(e, sphere);
    ms_hrtf_release(sphere);

    ms_physics_desc desc;
    ms_physics_desc_default(&desc);
    desc.voices = 4;
    desc.contacts = 4;
    ms_physics* physics = nullptr;
    check(ms_physics_create(&desc, &physics) == MS_OK, "physics created");
    if (!physics) {
        ++failures;
        ms_engine_destroy(e);
        return;
    }
    ms_physics_material steel;
    ms_physics_material_preset("steel", &steel);
    ms_modal_model plate = 0;
    ms_physics_model_box(physics, &steel, ms_vec3{0.6f, 0.4f, 0.004f}, 0, &plate);
    ms_engine_set_block_callback(e, ms_physics_block_callback, physics);
    ms_engine_commit(e);

    std::vector<float> left(kBlock), right(kBlock);
    ms_engine_render(e, left.data(), right.data(), kBlock);
    ms_physics_impact_desc hit{};
    hit.instance = 5;
    hit.model = plate;
    hit.impulse_ns = 0.1f;
    hit.contact_seconds = 1e-4f;
    hit.weight = 1.0f;
    check(ms_physics_impact(physics, &hit) == MS_OK, "impact posted");
    ms_engine_render(e, left.data(), right.data(), kBlock);
    check(ms_physics_voice_instance(physics, 0) == 5, "an impact starts with no source playing");

    ms_source sources[4];
    for (std::uint32_t i = 0; i < 4; ++i) {
        ms_source_create(e, &sources[i]);
        ms_source_params p;
        ms_source_params_default(&p);
        p.position = {static_cast<float>(i) - 1.5f, 0.0f, -2.0f};
        p.doppler = 0.0f;
        p.propagation_delay = 1;  // A lead on the first pull makes the pulls of different voices different sizes.
        ms_source_set_params(e, sources[i], &p);
        ms_source_set_callback(e, sources[i], ms_physics_source_callback, ms_physics_voice(physics, i));
        ms_source_play(e, sources[i], 0);
    }
    ms_engine_commit(e);

    ms_acoustics_config config;
    ms_acoustics_config_default(&config);
    ms_acoustics* acoustics = ms_acoustics_create(&config);
    const ms_acoustic_box bounds{{0, 0, 0}, {8, 6, 3}};
    ms_acoustic_room room{};
    room.bounds = &bounds;
    room.bound_count = 1;
    room.probe = {4, 3, 1.5f};
    room.volume = 144.0f;
    room.surface_area = 180.0f;
    for (float& area : room.absorption_area) area = 40.0f;
    ms_acoustics_room_add(acoustics, &room);
    check(ms_acoustics_reverb_prepare(acoustics) == MS_OK, "reverb prepared");
    std::vector<float> send(kBlock, 0.1f), wetL(kBlock), wetR(kBlock);
    const float* inputs[] = {send.data()};
    float* lefts[] = {wetL.data()};
    float* rights[] = {wetR.data()};

    double energy = 0.0;
    rt_probe_reset();
    rt_probe_watch_this_thread();
    for (int block = 0; block < 1500; ++block) {
        ms_engine_render(e, left.data(), right.data(), kBlock);
        ms_acoustics_reverb_process(acoustics, inputs, lefts, rights, kBlock);
        for (std::uint32_t i = 0; i < kBlock; ++i) energy += static_cast<double>(left[i]) * left[i];
    }
    rt_probe_stop();
    const unsigned long long forbidden = rt_probe_count();
    std::printf("  1500 blocks of stream pulls, physics clock and late reverb: %llu forbidden calls\n", forbidden);
    check(forbidden == 0, "the physics callback, stream pulls and reverb neither allocate nor block");
    check(energy > 0.0 && std::isfinite(energy), "the impact is heard through the engine");

    ms_engine_stats stats;
    ms_engine_get_stats(e, &stats);
    check(stats.faults == 0, "no faults");
    ms_engine_destroy(e);
    ms_physics_destroy(physics);
    ms_acoustics_destroy(acoustics);
}

}  // namespace

int main() {
    tapDelayIsTheGap();
    portalDelayIsTheGap();
    physicsOnTheBlockClock();
    std::printf(failures ? "marathon spatial integration: %d FAILED\n" : "marathon spatial integration: all passed\n", failures);
    return failures ? 1 : 0;
}
