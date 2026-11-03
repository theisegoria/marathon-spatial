// SPDX-License-Identifier: Apache-2.0
// Offline tests of the core. They link marathon_spatial_core alone, so nothing here can reach an audio device, and every
// signal goes into memory. Numbers print with their tolerance; a failure names what was measured.
#include <Accelerate/Accelerate.h>
#include <CoreFoundation/CoreFoundation.h>
#include <os/lock.h>

#include <algorithm>
#include <atomic>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <mutex>
#include <random>
#include <string>
#include <thread>
#include <vector>

#include "RtProbe.h"
#include "core/Ambisonics.hpp"
#include "core/Dsp.hpp"
#include "core/Engine.hpp"
#include "core/Crosstalk.hpp"
#include "core/DeviceChannels.hpp"
#include "core/Fft.hpp"
#include "core/Speakers.hpp"
#include "core/WorkgroupJoin.hpp"
#include "hrtf/HrtfSet.hpp"
#include "marathon_spatial/marathon_spatial.h"

using namespace mspatial;

namespace {

int gFailures = 0, gChecks = 0;
const char* gTest = "";

#define CHECK(cond, ...)                                                       \
    do {                                                                       \
        ++gChecks;                                                             \
        if (!(cond)) {                                                         \
            ++gFailures;                                                       \
            std::printf("  FAIL [%s] %s:%d %s: ", gTest, __FILE__, __LINE__, #cond); \
            std::printf(__VA_ARGS__);                                          \
            std::printf("\n");                                                 \
        }                                                                      \
    } while (0)
#define OK(expr) CHECK((expr) == MS_OK, "status %s", ms_status_string(expr))

constexpr float kRate = 48000.0f;
constexpr double kSpeed = 343.0;
constexpr float kRadius = 0.0875f;

struct Render {
    std::vector<float> l, r;
};

ms_engine* makeEngine(ms_hrtf* hrtf, std::uint32_t block = 128, std::uint32_t order = 3, float smoothing = 0.0f, std::uint32_t voices = 64) {
    ms_engine_desc d;
    ms_engine_desc_default(&d);
    d.block_frames = block;
    d.ambisonic_order = order;
    d.smoothing_seconds = smoothing;
    d.max_hrtf_voices = voices;
    ms_engine* e = nullptr;
    if (ms_engine_create(&d, &e) != MS_OK) std::abort();
    if (hrtf) ms_engine_set_hrtf(e, hrtf);
    return e;
}

Render renderFrames(ms_engine* e, std::uint32_t frames) {
    Render out;
    out.l.resize(frames);
    out.r.resize(frames);
    ms_engine_render(e, out.l.data(), out.r.data(), frames);
    return out;
}

// A plain distance-free configuration: what a head-model measurement wants.
ms_source_params plainParams(ms_vec3 position) {
    ms_source_params p;
    ms_source_params_default(&p);
    p.position = position;
    p.distance_model = MS_DISTANCE_NONE;
    p.air_absorption = 0.0f;
    p.near_field = 0.0f;
    p.doppler = 0.0f;
    return p;
}

ms_vec3 horizontal(double azimuthDegrees, float distance) {  // Positive azimuth is to the right.
    const double a = azimuthDegrees * M_PI / 180.0;
    return {static_cast<float>(std::sin(a)) * distance, 0.0f, static_cast<float>(-std::cos(a)) * distance};
}

Render impulseThrough(ms_hrtf* hrtf, ms_vec3 position, std::uint32_t frames = 2048) {
    ms_engine* e = makeEngine(hrtf);
    std::vector<float> impulse(64, 0.0f);
    impulse[0] = 1.0f;
    ms_buffer b;
    ms_source s;
    OK(ms_buffer_create(e, impulse.data(), impulse.size(), 1, kRate, &b));
    OK(ms_source_create(e, &s));
    const ms_source_params p = plainParams(position);
    OK(ms_source_set_params(e, s, &p));
    OK(ms_source_set_buffer(e, s, b, 0, 0));
    OK(ms_source_play(e, s, 0));
    OK(ms_engine_commit(e));
    Render r = renderFrames(e, frames);
    ms_engine_destroy(e);
    return r;
}

// Magnitude response of an impulse response at a frequency, by direct DFT.
double magnitudeDb(const std::vector<float>& x, double hz) {
    std::complex<double> sum = 0.0;
    for (std::size_t n = 0; n < x.size(); ++n) sum += static_cast<double>(x[n]) * std::polar(1.0, -2.0 * M_PI * hz * static_cast<double>(n) / kRate);
    return 20.0 * std::log10(std::abs(sum) + 1e-20);
}

// Average over a band, in power.
double bandDb(const std::vector<float>& x, double low, double high) {
    double total = 0.0;
    constexpr int kSteps = 24;
    for (int i = 0; i < kSteps; ++i) {
        const double hz = low * std::pow(high / low, (i + 0.5) / kSteps);
        total += std::pow(10.0, magnitudeDb(x, hz) / 10.0);
    }
    return 10.0 * std::log10(total / kSteps);
}

// Time of a band's arrival, seconds: the peak of the envelope's interaural cross-correlation.
double interauralDelaySeconds(const std::vector<float>& l, const std::vector<float>& r, double cutoffHz) {
    const auto lowpass = [&](const std::vector<float>& in) {
        std::vector<float> out = in;
        core::Biquad a = core::designLowpass(static_cast<float>(cutoffHz), 0.7071f, kRate), b = a;
        a.process(out.data(), static_cast<std::uint32_t>(out.size()));
        b.process(out.data(), static_cast<std::uint32_t>(out.size()));
        return out;
    };
    const std::vector<float> fl = lowpass(l), fr = lowpass(r);
    constexpr int kReach = 64;
    const auto n = static_cast<int>(fl.size());
    std::vector<double> correlation(2 * kReach + 1);
    for (int lag = -kReach; lag <= kReach; ++lag) {
        double sum = 0.0;
        for (int i = std::max(0, -lag); i < std::min(n, n - lag); ++i) sum += static_cast<double>(fl[static_cast<std::size_t>(i)]) * fr[static_cast<std::size_t>(i + lag)];
        correlation[static_cast<std::size_t>(lag + kReach)] = sum;
    }
    const auto best = std::distance(correlation.begin(), std::max_element(correlation.begin(), correlation.end()));
    double offset = 0.0;
    if (best > 0 && best < 2 * kReach) {
        const double a = correlation[static_cast<std::size_t>(best - 1)], b = correlation[static_cast<std::size_t>(best)], c = correlation[static_cast<std::size_t>(best + 1)];
        offset = 0.5 * (a - c) / (a - 2.0 * b + c);
    }
    return (static_cast<double>(best - kReach) + offset) / kRate;  // Positive: right ear later.
}

// First crossing of 10 percent of the channel's peak, with linear interpolation, in seconds. The ear's onset.
double onsetSeconds(const std::vector<float>& x) {
    float peak = 0.0f;
    for (float v : x) peak = std::max(peak, std::fabs(v));
    const float threshold = 0.1f * peak;
    for (std::size_t i = 1; i < x.size(); ++i)
        if (std::fabs(x[i]) >= threshold) {
            const double a = std::fabs(x[i - 1]), b = std::fabs(x[i]);
            return (static_cast<double>(i - 1) + (threshold - a) / (b - a + 1e-20)) / kRate;
        }
    return 0.0;
}

double sphereDb(double rho, double cosAngle, double hz, double radius = kRadius) {
    double re, im;
    hrtf::sphereTransfer(rho, cosAngle, 2.0 * M_PI * hz * radius / kSpeed, re, im);
    return 10.0 * std::log10(re * re + im * im);
}

std::string dataPath(const char* name) { return std::string(MS_DATA_DIR) + "/" + name; }

ms_hrtf* gSphere = nullptr;
std::vector<std::pair<std::string, ms_hrtf*>> gMeasured;

void testSphereModel() {
    gTest = "sphere model";
    // Far away and at low frequency the sphere is acoustically transparent.
    CHECK(std::fabs(sphereDb(2000.0, 1.0, 20.0)) < 0.05, "low frequency %f dB", sphereDb(2000.0, 1.0, 20.0));
    CHECK(std::fabs(sphereDb(2000.0, -1.0, 20.0)) < 0.05, "low frequency rear");
    // The static limit has a closed form: Sum (2n + 1) / (n + 1) rho^-n on the near pole.
    double staticLimit = 0.0;
    for (int n = 0; n < 200; ++n) staticLimit += (2.0 * n + 1.0) / (n + 1.0) * std::pow(1.0 / 2.0, n);
    CHECK(std::fabs(sphereDb(2.0, 1.0, 1.0) - 20.0 * std::log10(staticLimit)) < 0.02, "static limit %f dB, closed form %f", sphereDb(2.0, 1.0, 1.0), 20.0 * std::log10(staticLimit));
    // Pressure doubling on the lit pole as ka grows.
    const double lit = sphereDb(200.0, 1.0, 16000.0);
    CHECK(std::fabs(lit - 6.02) < 1.0, "lit pole at high ka %f dB", lit);
    // The bright spot: at the antipode, in phase round both ways, more than 10 dB above the flanks at ka near 10.
    const double spot = sphereDb(200.0, -1.0, 6000.0), flank = sphereDb(200.0, std::cos(150.0 * M_PI / 180.0), 6000.0);
    CHECK(spot > flank + 8.0, "bright spot %f vs flank %f", spot, flank);
    // A shadowed point falls with frequency, the classic head shadow.
    double deepest = 0.0;
    for (int angle = 90; angle < 180; ++angle) deepest = std::min(deepest, sphereDb(200.0, std::cos(angle * M_PI / 180.0), 8000.0));
    CHECK(deepest < -8.0, "shadow %f", deepest);
    // Phase: delay against centre arrival is a / c early at the lit pole, at high frequency.
    double re1, im1, re2, im2;
    const double hz = 12000.0, d = 50.0;
    hrtf::sphereTransfer(200.0, 1.0, 2 * M_PI * hz * kRadius / kSpeed, re1, im1);
    hrtf::sphereTransfer(200.0, 1.0, 2 * M_PI * (hz + d) * kRadius / kSpeed, re2, im2);
    const double phase = std::arg(std::complex<double>(re2, im2) * std::conj(std::complex<double>(re1, im1)));
    const double groupDelay = -phase / (2 * M_PI * d);
    CHECK(std::fabs(groupDelay + kRadius / kSpeed) < 0.15 * kRadius / kSpeed, "group delay %f us, expected %f", groupDelay * 1e6, -kRadius / kSpeed * 1e6);
}

void testFilterIdentityAndLevels() {
    gTest = "levels";
    // The engine's convolution scale: the rendered response must equal the set's stored filter, at unit gain.
    const Render front = impulseThrough(gSphere, horizontal(0.0, 1.5));
    double energy = 0.0;
    for (std::size_t i = 0; i < front.l.size(); ++i) energy += static_cast<double>(front.l[i]) * front.l[i];
    const double db = magnitudeDb(front.l, 100.0);
    CHECK(std::fabs(db) < 0.5, "100 Hz gain %f dB", db);
    std::printf("  front ear gain at 100 Hz %.3f dB, energy %.4f\n", db, energy);
}

void testHeadModelCues(const char* label, ms_hrtf* hrtf, bool analytic) {
    gTest = "head model cues";
    ms_hrtf_info info;
    ms_hrtf_get_info(hrtf, &info);
    const double a = info.head_radius_m;
    std::printf("  %s: %u directions, head radius %.1f mm\n", label, info.directions, a * 1000.0);
    double worstItd = 0.0, previous = -1.0;
    std::printf("    az   onset ITD us  Woodworth   xcorr<1k us  ILD 1-2k  ILD 3-6k  ILD 8-12k\n");
    for (int azimuth = 0; azimuth <= 90; azimuth += 15) {
        const Render r = impulseThrough(hrtf, horizontal(azimuth, 1.5));
        const double onset = (onsetSeconds(r.l) - onsetSeconds(r.r)) * 1e6;
        const double model = hrtf::woodworthItdSeconds(azimuth * M_PI / 180.0, a, kSpeed) * 1e6;
        const double correlated = -interauralDelaySeconds(r.r, r.l, 1000.0) * 1e6;
        const double ild1 = bandDb(r.r, 1000, 2000) - bandDb(r.l, 1000, 2000);
        const double ild2 = bandDb(r.r, 3000, 6000) - bandDb(r.l, 3000, 6000);
        const double ild3 = bandDb(r.r, 8000, 12000) - bandDb(r.l, 8000, 12000);
        std::printf("    %3d  %10.1f  %9.1f  %11.1f  %8.2f  %8.2f  %9.2f\n", azimuth, onset, model, correlated, ild1, ild2, ild3);
        worstItd = std::max(worstItd, std::fabs(onset - model));
        CHECK(onset >= previous - 25.0, "ITD not monotonic at %d: %f after %f", azimuth, onset, previous);
        previous = onset;
        if (analytic) {
            const double angleRight = std::cos((90.0 - azimuth) * M_PI / 180.0), angleLeft = std::cos((90.0 + azimuth) * M_PI / 180.0);
            const double rho = info.measured_distance_m / a;
            for (double hz : {1000.0, 2000.0, 4000.0, 8000.0}) {
                const double modelIld = sphereDb(rho, angleRight, hz, a) - sphereDb(rho, angleLeft, hz, a);
                const double measured = magnitudeDb(r.r, hz) - magnitudeDb(r.l, hz);
                CHECK(std::fabs(modelIld - measured) < 1.5, "ILD at %d deg %.0f Hz: %f, sphere %f", azimuth, hz, measured, modelIld);
            }
        }
    }
    std::printf("    worst onset ITD error against Woodworth: %.1f us\n", worstItd);
    if (analytic) CHECK(worstItd < 45.0, "worst %f us", worstItd);
    else CHECK(worstItd < 130.0, "worst %f us", worstItd);
}

double mirrorError(const std::vector<float>& a, const std::vector<float>& b) {
    double difference = 0.0, reference = 1e-30;
    for (std::size_t i = 0; i < a.size(); ++i) {
        const double d = static_cast<double>(a[i]) - b[i];
        difference += d * d;
        reference += static_cast<double>(a[i]) * a[i];
    }
    return 10.0 * std::log10(difference / reference + 1e-30);
}

void testSymmetry() {
    gTest = "impulse symmetry";
    double worst = -300.0;
    for (int azimuth : {15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165}) {
        const Render plus = impulseThrough(gSphere, horizontal(azimuth, 1.5)), minus = impulseThrough(gSphere, horizontal(-azimuth, 1.5));
        worst = std::max({worst, mirrorError(plus.l, minus.r), mirrorError(plus.r, minus.l)});
    }
    std::printf("  worst mirror error, sphere set: %.1f dB\n", worst);
    CHECK(worst < -80.0, "mirror error %f dB", worst);
    for (const auto& [name, set] : gMeasured) {
        const Render plus = impulseThrough(set, horizontal(40, 1.5)), minus = impulseThrough(set, horizontal(-40, 1.5));
        std::printf("  %s mirror difference at 40 degrees, a property of the measured head: %.1f dB\n", name.c_str(), mirrorError(plus.l, minus.r));
    }
    // Turning the listener and placing the source are the same scene: 30 degrees right of forward either way.
    std::vector<float> noise(6000);
    std::mt19937 rng(7);
    for (float& v : noise) v = std::uniform_real_distribution<float>(-1, 1)(rng);
    Render results[2];
    for (int variant = 0; variant < 2; ++variant) {
        ms_engine* e = makeEngine(gSphere);
        ms_buffer buffer;
        ms_source s;
        OK(ms_buffer_create(e, noise.data(), noise.size(), 1, kRate, &buffer));
        OK(ms_source_create(e, &s));
        const ms_vec3 offset = {3.0f, -1.0f, 2.0f};
        ms_vec3 at = variant == 0 ? horizontal(30.0, 2.0) : horizontal(70.0, 2.0);
        if (variant == 1) {
            at = {at.x + offset.x, at.y + offset.y, at.z + offset.z};
            OK(ms_listener_set_pose(e, offset, {0, 0, 0}, horizontal(40.0, 1.0), {0, 1, 0}));
        }
        const ms_source_params p = plainParams(at);
        OK(ms_source_set_params(e, s, &p));
        OK(ms_source_set_buffer(e, s, buffer, 0, 0));
        OK(ms_source_play(e, s, 0));
        OK(ms_engine_commit(e));
        results[variant] = renderFrames(e, 6000);
        ms_engine_destroy(e);
    }
    const double rotated = mirrorError(results[0].l, results[1].l);
    std::printf("  listener turned against source moved: %.1f dB\n", rotated);
    CHECK(rotated < -60.0, "%f dB", rotated);
}

// Everything at once, for the allocation and determinism tests: moving sources through both paths, a callback, a bed,
// reflections, buses, a send and a processor.
struct Scene {
    ms_engine* e{};
    ms_source sources[24]{};
    ms_bus bus{};
    std::uint32_t phase{};
    static uint32_t tone(void* user, float* mono, uint32_t frames) {
        auto* phase = static_cast<std::uint32_t*>(user);
        for (uint32_t i = 0; i < frames; ++i) mono[i] = 0.2f * std::sin(0.03f * static_cast<float>((*phase)++));
        return frames;
    }
    static void processor(void*, float* left, float* right, uint32_t frames, uint64_t) {
        for (uint32_t i = 0; i < frames; ++i) {
            left[i] *= 0.9f;
            right[i] *= 0.9f;
        }
    }
    void build(ms_hrtf* hrtf, std::uint32_t voices, std::uint32_t block = 128) {
        e = makeEngine(hrtf, block, 3, 0.005f, voices);
        std::vector<float> noise(24000), bedData(16 * 24000);
        std::mt19937 rng(11);
        for (float& v : noise) v = std::uniform_real_distribution<float>(-0.3f, 0.3f)(rng);
        for (float& v : bedData) v = std::uniform_real_distribution<float>(-0.05f, 0.05f)(rng);
        ms_buffer buffer, bedBuffer;
        OK(ms_buffer_create(e, noise.data(), noise.size(), 1, kRate, &buffer));
        OK(ms_buffer_create(e, bedData.data(), 24000, 16, kRate, &bedBuffer));
        OK(ms_bus_create(e, MS_BUS_MASTER, &bus));
        OK(ms_bus_set_gain(e, bus, 0.8f));
        OK(ms_bus_set_filter(e, bus, 40.0f, 15000.0f));
        OK(ms_bus_set_processor(e, bus, processor, nullptr));
        for (int i = 0; i < 24; ++i) {
            OK(ms_source_create(e, &sources[i]));
            ms_source_params p;
            ms_source_params_default(&p);
            p.bus = bus;
            p.send_bus[0] = i % 2 ? MS_BUS_MASTER : 0;
            p.send_gain[0] = 0.1f;
            p.directivity_alpha = 0.5f;
            p.directivity_hf_db = -6.0f;
            p.lowpass_hz = 9000.0f;
            p.eq_high_db = -3.0f;
            p.position = {static_cast<float>(i) * 0.5f - 6.0f, 0.0f, -2.0f};
            OK(ms_source_set_params(e, sources[i], &p));
            if (i == 0) OK(ms_source_set_callback(e, sources[i], tone, &phase));
            else OK(ms_source_set_buffer(e, sources[i], buffer, 0, 1));
            OK(ms_source_play(e, sources[i], 0));
        }
        ms_reflection_tap taps[6];
        for (int t = 0; t < 6; ++t) taps[t] = {0.01f + 0.01f * static_cast<float>(t), {1.0f, 0.0f, static_cast<float>(t) - 3.0f}, {0.5f, 0.3f, 0.2f}};
        for (int i = 0; i < 4; ++i) OK(ms_source_set_reflections(e, sources[i], taps, 6));
        ms_bed bed;
        OK(ms_bed_create(e, &bed));
        OK(ms_bed_set_buffer(e, bed, bedBuffer, 1));
        OK(ms_bed_play(e, bed, 0));
        OK(ms_engine_commit(e));
    }
    void move(int step) {
        const float t = static_cast<float>(step) * 0.05f;
        for (int i = 0; i < 24; ++i) {
            const float a = t + static_cast<float>(i) * 0.3f, radius = 0.5f + 0.3f * static_cast<float>(i);
            OK(ms_source_set_pose(e, sources[i], {radius * std::sin(a), 0.3f * std::sin(t), -radius * std::cos(a)}, {radius * 3.0f * std::cos(a), 0.0f, radius * 3.0f * std::sin(a)}));
        }
        OK(ms_listener_set_pose(e, {0.2f * std::sin(t), 0, 0}, {0.1f, 0, 0}, horizontal(10.0 * std::sin(t), 1.0), {0, 1, 0}));
        const float half = 0.3f * std::sin(t * 1.3f);
        ms_head_set_orientation(e, {0, std::sin(half), 0, std::cos(half)});
        OK(ms_engine_commit(e));
    }
};

void testNoAllocationOnRenderSide() {
    gTest = "no allocation";
    // The probe first has to prove it can see what it is looking for.
    rt_probe_reset();
    rt_probe_watch_this_thread();
    const auto expectSeen = [](const char* what, const std::function<void()>& action) {
        const unsigned long long before = rt_probe_count();
        action();
        const unsigned long long after = rt_probe_count();
        std::printf("  probe control, %s: %llu\n", what, after - before);
#ifndef MS_RT_PROBE_OFF
        CHECK(after > before, "probe blind to %s", what);
#endif
    };
    expectSeen("malloc", [] {
        void* volatile p = std::malloc(32);
        std::free(p);
    });
    expectSeen("operator new", [] {
        std::vector<int>* volatile v = new std::vector<int>(100);
        delete v;
    });
    expectSeen("an allocation inside a system library", [] {
        CFStringRef text = CFStringCreateWithCString(nullptr, "a string long enough to need its own storage", kCFStringEncodingUTF8);
        CFRelease(text);
    });
    expectSeen("std::mutex", [] {
        std::mutex mutex;
        mutex.lock();
        mutex.unlock();
    });
    expectSeen("os_unfair_lock", [] {
        os_unfair_lock lock = OS_UNFAIR_LOCK_INIT;
        os_unfair_lock_lock(&lock);
        os_unfair_lock_unlock(&lock);
    });
    rt_probe_stop();

    Scene scene;
    scene.build(gSphere, 8);
    std::vector<float> l(128), r(128);
    rt_probe_reset();
    std::atomic<bool> running{true};
    std::atomic<std::uint64_t> blocks{0};
    std::thread audio([&] {
        rt_probe_watch_this_thread();
        while (running.load(std::memory_order_relaxed)) {
            ms_engine_render(scene.e, l.data(), r.data(), 128);
            blocks.fetch_add(1, std::memory_order_relaxed);
        }
        rt_probe_stop();
    });
    for (int step = 0; step < 400; ++step) {
        scene.move(step);
        if (step % 40 == 0) {  // Structural churn while the render side runs.
            OK(ms_source_stop(scene.e, scene.sources[step / 40 + 2]));
            OK(ms_source_play(scene.e, scene.sources[(step / 40 + 1) % 24 + 1 > 23 ? 2 : step / 40 + 1], 0));
            OK(ms_engine_commit(scene.e));
        }
        std::this_thread::sleep_for(std::chrono::microseconds(200));
    }
    while (blocks.load(std::memory_order_relaxed) < 200) std::this_thread::yield();
    running = false;
    audio.join();
    std::printf("  %llu blocks with commits racing: %llu forbidden calls on the render thread\n", static_cast<unsigned long long>(blocks.load()), rt_probe_count());
    CHECK(rt_probe_count() == 0, "%llu forbidden calls", rt_probe_count());
    CHECK(blocks.load() > 100, "only %llu blocks", static_cast<unsigned long long>(blocks.load()));
    ms_engine_destroy(scene.e);
}

std::uint64_t hashOf(const Render& r) {
    std::uint64_t h = 14695981039346656037ull;
    for (const auto* v : {&r.l, &r.r}) {
        const auto* bytes = reinterpret_cast<const unsigned char*>(v->data());
        for (std::size_t i = 0; i < v->size() * sizeof(float); ++i) h = (h ^ bytes[i]) * 1099511628211ull;
    }
    return h;
}

Render deterministicRun(const std::vector<std::uint32_t>& chunks) {
    Scene scene;
    scene.build(gSphere, 8);
    Render out;
    constexpr std::uint32_t kFrames = 48000;
    out.l.resize(kFrames);
    out.r.resize(kFrames);
    std::uint32_t done = 0, chunkAt = 0;
    int step = 0;
    std::uint32_t nextMove = 0;
    while (done < kFrames) {
        // Control happens at frame positions, independent of chunking.
        if (done >= nextMove) {
            scene.move(step++);
            nextMove += 1024;
        }
        const std::uint32_t untilMove = nextMove - done;
        const std::uint32_t chunk = std::min({chunks[chunkAt++ % chunks.size()], kFrames - done, untilMove});
        ms_engine_render(scene.e, out.l.data() + done, out.r.data() + done, chunk);
        done += chunk;
    }
    ms_engine_destroy(scene.e);
    return out;
}

void testDeterminism() {
    gTest = "determinism";
    const Render first = deterministicRun({128}), second = deterministicRun({128});
    CHECK(hashOf(first) == hashOf(second), "two runs differ");
    // Chunking changes when a commit is seen unless blocks align, so the moves above sit on block boundaries for these sizes.
    const Render chunked = deterministicRun({64, 64, 256});
    const Render uneven = deterministicRun({128, 1024, 256, 512});
    CHECK(hashOf(first) == hashOf(chunked), "chunked differs");
    CHECK(hashOf(first) == hashOf(uneven), "uneven differs");
    std::printf("  hash %016llx over one second of a busy scene, four runs identical in every bit\n", static_cast<unsigned long long>(hashOf(first)));
    double energy = 0.0;
    for (float v : first.l) energy += static_cast<double>(v) * v;
    CHECK(energy > 1.0, "scene is silent: energy %f", energy);
}

double dominantFrequency(const std::vector<float>& x, std::size_t begin, std::size_t end) {  // Mean period by interpolated rising zero crossings.
    double first = -1.0, last = -1.0;
    int crossings = 0;
    for (std::size_t i = begin + 1; i < end; ++i)
        if (x[i - 1] < 0.0f && x[i] >= 0.0f) {
            const double t = static_cast<double>(i - 1) + static_cast<double>(-x[i - 1]) / (static_cast<double>(x[i]) - x[i - 1]);
            if (first < 0.0) first = t;
            last = t;
            ++crossings;
        }
    return crossings > 1 ? (crossings - 1) * kRate / (last - first) : 0.0;
}

void testDoppler() {
    gTest = "doppler";
    std::vector<float> tone(96000);
    for (std::size_t i = 0; i < tone.size(); ++i) tone[i] = 0.5f * std::sin(2.0f * core::kPi * 1000.0f * static_cast<float>(i) / kRate);
    for (const float speed : {-34.3f, 34.3f}) {
        ms_engine* e = makeEngine(gSphere);
        ms_buffer b;
        ms_source s;
        OK(ms_buffer_create(e, tone.data(), tone.size(), 1, kRate, &b));
        OK(ms_source_create(e, &s));
        ms_source_params p = plainParams({0, 0, -10});
        p.doppler = 1.0f;
        p.velocity = {0, 0, -speed};  // Positive speed recedes from the listener at the origin, along -z.
        OK(ms_source_set_params(e, s, &p));
        OK(ms_source_set_buffer(e, s, b, 0, 1));
        OK(ms_source_play(e, s, 0));
        OK(ms_engine_commit(e));
        const Render r = renderFrames(e, 48000);
        const double measured = dominantFrequency(r.l, 4800, 48000);
        const double expected = 1000.0 * kSpeed / (kSpeed + speed);
        std::printf("  source at %+.1f m/s along the line: %.2f Hz, expected %.2f\n", speed, measured, expected);
        CHECK(std::fabs(measured - expected) < 0.5, "%f vs %f", measured, expected);
        ms_engine_destroy(e);
    }
}

void testDistanceAndAir() {
    gTest = "distance and air";
    ms_engine* e = makeEngine(gSphere);
    std::vector<float> impulse(32, 0.0f);
    impulse[0] = 1.0f;
    ms_engine_destroy(e);
    const auto response = [&](float distance, float air, std::uint32_t model) {
        ms_engine* engine = makeEngine(gSphere);
        ms_buffer b;
        ms_source s;
        OK(ms_buffer_create(engine, impulse.data(), impulse.size(), 1, kRate, &b));
        OK(ms_source_create(engine, &s));
        ms_source_params p = plainParams({0, 0, -distance});
        p.air_absorption = air;
        p.distance_model = model;
        p.ref_distance = 1.0f;
        p.max_distance = 10000.0f;
        OK(ms_source_set_params(engine, s, &p));
        OK(ms_source_set_buffer(engine, s, b, 0, 0));
        OK(ms_source_play(engine, s, 0));
        OK(ms_engine_commit(engine));
        Render r = renderFrames(engine, 4096);
        ms_engine_destroy(engine);
        return r;
    };
    const Render near = response(1.0f, 0.0f, MS_DISTANCE_INVERSE), far = response(8.0f, 0.0f, MS_DISTANCE_INVERSE);
    const double law = magnitudeDb(near.l, 500.0) - magnitudeDb(far.l, 500.0);
    CHECK(std::fabs(law - 18.06) < 0.1, "inverse law %f dB over 3 doublings", law);
    const Render dry = response(100.0f, 0.0f, MS_DISTANCE_NONE), wet = response(100.0f, 1.0f, MS_DISTANCE_NONE);
    std::printf("  air over 100 m at 20 C, 50%%RH:   Hz  rendered  ISO 9613-1\n");
    double worst = 0.0;
    for (double hz : {1000.0, 2000.0, 4000.0, 8000.0, 12000.0, 16000.0}) {
        const double rendered = magnitudeDb(dry.l, hz) - magnitudeDb(wet.l, hz);
        const double iso = 100.0 * core::airAttenuationDbPerMetre(hz, 20.0, 50.0);
        std::printf("    %6.0f  %7.2f  %7.2f\n", hz, rendered, iso);
        worst = std::max(worst, std::fabs(rendered - iso));
    }
    std::printf("  worst air error over 100 m: %.2f dB\n", worst);
    CHECK(worst < 6.0, "%f dB", worst);
    // ISO 9613-1 table values at 20 C, 50 %RH, 1 atm: 4.66 dB/km at 1 kHz, 22.7 at 4 kHz (rounded in the standard's table).
    CHECK(std::fabs(core::airAttenuationDbPerMetre(1000.0, 20.0, 50.0) * 1000.0 - 4.66) < 0.3, "ISO at 1 kHz %f", core::airAttenuationDbPerMetre(1000.0, 20.0, 50.0) * 1000.0);
}

void testNearField() {
    gTest = "near field";
    const float rho = 2.5f, distance = kRadius * rho;
    ms_engine* engine = makeEngine(gSphere);
    ms_engine_destroy(engine);
    std::vector<float> impulse(32, 0.0f);
    impulse[0] = 1.0f;
    const auto response = [&](float nearField) {
        ms_engine* e = makeEngine(gSphere);
        ms_buffer b;
        ms_source s;
        OK(ms_buffer_create(e, impulse.data(), impulse.size(), 1, kRate, &b));
        OK(ms_source_create(e, &s));
        ms_source_params p = plainParams(horizontal(90.0, distance));
        p.near_field = nearField;
        OK(ms_source_set_params(e, s, &p));
        OK(ms_source_set_buffer(e, s, b, 0, 0));
        OK(ms_source_play(e, s, 0));
        OK(ms_engine_commit(e));
        Render r = renderFrames(e, 4096);
        ms_engine_destroy(e);
        return r;
    };
    const Render on = response(1.0f), off = response(0.0f);
    const double rho0 = 1.5 / kRadius;
    std::printf("  source at %.1f cm, 90 degrees:   Hz  ILD rendered  ILD sphere  ILD far set\n", distance * 100.0);
    double worst = 0.0;
    for (double hz : {200.0, 500.0, 1000.0, 2000.0}) {
        const double rendered = magnitudeDb(on.r, hz) - magnitudeDb(on.l, hz);
        const double model = sphereDb(rho, 1.0, hz) - sphereDb(rho, -1.0, hz);
        const double far = magnitudeDb(off.r, hz) - magnitudeDb(off.l, hz);
        std::printf("    %6.0f  %8.2f  %8.2f  %8.2f\n", hz, rendered, model, far);
        CHECK(rendered > far + 1.5, "no near-field rise at %f Hz", hz);
        worst = std::max(worst, std::fabs(rendered - model));
    }
    (void)rho0;
    CHECK(worst < 3.5, "worst %f dB", worst);
}

void testAmbisonics() {
    gTest = "ambisonics";
    float gains[core::kMaxChannels], rotated[core::kMaxChannels], expected[core::kMaxChannels];
    // SN3D and ACN by their definitions.
    core::encodeDirection({0.6f, 0.0f, 0.8f}, 2, gains);
    CHECK(std::fabs(gains[0] - 1.0f) < 1e-6f, "W");
    CHECK(std::fabs(gains[1]) < 1e-6f && std::fabs(gains[2] - 0.8f) < 1e-6f && std::fabs(gains[3] - 0.6f) < 1e-6f, "Y Z X");
    CHECK(std::fabs(gains[6] - 0.5f * (3.0f * 0.64f - 1.0f)) < 1e-6f, "R");
    CHECK(std::fabs(gains[7] - std::sqrt(3.0f) * 0.8f * 0.6f) < 1e-6f, "S");
    // Rotation of the channels equals encoding the rotated direction, to order 7.
    core::ShRotation rotation;
    std::mt19937 rng(3);
    std::normal_distribution<float> normal;
    double worst = 0.0;
    std::vector<float> matrix(core::kRotationCells);
    for (int trial = 0; trial < 200; ++trial) {
        const core::Quat q = core::normalized(core::Quat{normal(rng), normal(rng), normal(rng), normal(rng)});
        const core::Vec3 axisX = core::rotate(q, {1, 0, 0}), axisY = core::rotate(q, {0, 1, 0}), axisZ = core::rotate(q, {0, 0, 1});
        const float m[9] = {axisX.x, axisY.x, axisZ.x, axisX.y, axisY.y, axisZ.y, axisX.z, axisY.z, axisZ.z};
        rotation.compute(m, 7, matrix.data());
        const core::Vec3 d = core::normalized(core::Vec3{normal(rng), normal(rng), normal(rng)});
        core::encodeDirection(d, 7, gains);
        core::encodeDirection(core::rotate(q, d), 7, expected);
        std::size_t cell = 0;
        for (std::uint32_t l = 0; l <= 7; ++l) {
            const std::uint32_t size = 2 * l + 1;
            for (std::uint32_t i = 0; i < size; ++i) {
                float sum = 0.0f;
                for (std::uint32_t j = 0; j < size; ++j, ++cell) sum += matrix[cell] * gains[l * l + j];
                rotated[l * l + i] = sum;
            }
        }
        for (std::uint32_t c = 0; c < 64; ++c) worst = std::max(worst, static_cast<double>(std::fabs(rotated[c] - expected[c])));
    }
    std::printf("  worst channel error of the rotation over 200 random cases, order 7: %.2e\n", worst);
    CHECK(worst < 2e-4, "%e", worst);

    // A plane wave in a bed against the same direction rendered directly, through the same set.
    constexpr std::uint32_t kFrames = 8192;
    for (std::uint32_t order : {1u, 3u, 5u}) {
        double worstIld = 0.0, worstItd = 0.0;
        for (int azimuth : {-60, 30, 90, 135}) {
            ms_engine* e = makeEngine(gSphere, 128, order);
            const std::uint32_t channels = (order + 1) * (order + 1);
            std::vector<float> data(static_cast<std::size_t>(kFrames) * channels, 0.0f);
            const ms_vec3 v = horizontal(azimuth, 1.0);
            std::vector<float> encode(channels);
            core::encodeDirection(core::headToAmbisonic(core::Vec3(v)), order, encode.data());
            std::vector<float> noise(32, 0.0f);  // An impulse: a flat spectrum, so band levels are the filters' own.
            noise[0] = 1.0f;
            for (std::uint32_t c = 0; c < channels; ++c) data[c] = encode[c];
            ms_buffer b;
            ms_bed bed;
            OK(ms_buffer_create(e, data.data(), kFrames, channels, kRate, &b));
            OK(ms_bed_create(e, &bed));
            OK(ms_bed_set_buffer(e, bed, b, 0));
            OK(ms_bed_play(e, bed, 0));
            OK(ms_engine_commit(e));
            const Render fromBed = renderFrames(e, kFrames);
            ms_engine_destroy(e);

            ms_engine* direct = makeEngine(gSphere, 128, 0);
            ms_buffer b2;
            ms_source s;
            OK(ms_buffer_create(direct, noise.data(), noise.size(), 1, kRate, &b2));
            OK(ms_source_create(direct, &s));
            const ms_source_params p = plainParams(horizontal(azimuth, 1.5));
            OK(ms_source_set_params(direct, s, &p));
            OK(ms_source_set_buffer(direct, s, b2, 0, 0));
            OK(ms_source_play(direct, s, 0));
            OK(ms_engine_commit(direct));
            const Render fromSource = renderFrames(direct, kFrames);
            ms_engine_destroy(direct);

            const auto ild = [](const Render& r, double a, double b) { return bandDb(r.r, a, b) - bandDb(r.l, a, b); };
            worstIld = std::max({worstIld, std::fabs(ild(fromBed, 500, 2000) - ild(fromSource, 500, 2000)), std::fabs(ild(fromBed, 3000, 8000) - ild(fromSource, 3000, 8000))});
            const double bedItd = interauralDelaySeconds(fromBed.l, fromBed.r, 1000.0) * 1e6, directItd = interauralDelaySeconds(fromSource.l, fromSource.r, 1000.0) * 1e6;
            worstItd = std::max(worstItd, std::fabs(bedItd - directItd));
            std::printf("    order %u az %4d: ILD 0.5-2k bed %6.2f direct %6.2f, 3-8k bed %6.2f direct %6.2f, ITD bed %7.1f direct %7.1f us\n", order, azimuth, ild(fromBed, 500, 2000), ild(fromSource, 500, 2000), ild(fromBed, 3000, 8000), ild(fromSource, 3000, 8000), bedItd, directItd);
        }
        std::printf("  order %u bed against direct: worst ILD difference %.2f dB, worst ITD difference under 1 kHz %.1f us\n", order, worstIld, worstItd);
        // Limits from what a truncated order can hold: first order cannot carry ITD past about 600 Hz, and the rigid sphere's
        // bright spot at exactly 90 degrees is a sharp feature that takes order 5 to resolve.
        CHECK(worstIld < (order == 1 ? 9.0 : order == 3 ? 5.0 : 2.5), "order %u ILD %f", order, worstIld);
        if (order > 1) CHECK(worstItd < (order == 3 ? 15.0 : 5.0), "order %u ITD %f", order, worstItd);
    }

    // World locked: turning the head turns the bed against it, by the same amount as turning the bed by hand the other way.
    Render turned[2];
    for (int variant = 0; variant < 2; ++variant) {
        ms_engine* e = makeEngine(gSphere, 128, 3);
        std::vector<float> data(static_cast<std::size_t>(4096) * 16);
        std::mt19937 local(9);
        for (float& x : data) x = std::uniform_real_distribution<float>(-0.1f, 0.1f)(local);
        ms_buffer b;
        ms_bed bed;
        OK(ms_buffer_create(e, data.data(), 4096, 16, kRate, &b));
        OK(ms_bed_create(e, &bed));
        OK(ms_bed_set_buffer(e, bed, b, 0));
        const float half = 0.4f;
        if (variant == 0) ms_head_set_orientation(e, {0, std::sin(half), 0, std::cos(half)});
        else OK(ms_bed_set_orientation(e, bed, {0, std::sin(-half), 0, std::cos(half)}, 0));
        OK(ms_bed_play(e, bed, 0));
        OK(ms_engine_commit(e));
        turned[variant] = renderFrames(e, 4096);
        ms_engine_destroy(e);
    }
    const double difference = mirrorError(turned[0].l, turned[1].l);
    std::printf("  head turned against bed turned: %.1f dB\n", difference);
    CHECK(difference < -60.0, "%f dB", difference);
}

// A sine sweeping around the head must keep its energy at the tone: the interpolation should not add clicks or zipper noise.
void testMovingSourceStaysClean() {
    gTest = "moving source";
    std::vector<float> tone(96000);
    for (std::size_t i = 0; i < tone.size(); ++i) tone[i] = 0.5f * std::sin(2.0f * core::kPi * 500.0f * static_cast<float>(i) / kRate);
    ms_engine* e = makeEngine(gSphere, 128, 3, 0.005f);
    ms_buffer b;
    ms_source s;
    OK(ms_buffer_create(e, tone.data(), tone.size(), 1, kRate, &b));
    OK(ms_source_create(e, &s));
    ms_source_params p = plainParams(horizontal(0.0, 2.0));
    OK(ms_source_set_params(e, s, &p));
    OK(ms_source_set_buffer(e, s, b, 0, 1));
    OK(ms_source_play(e, s, 0));
    OK(ms_engine_commit(e));
    constexpr std::uint32_t kBlocks = 3000, kHop = 128 * 6;  // A 60 Hz game feeding a circle every 3.3 s
    std::vector<float> left;
    left.reserve(kBlocks * 128);
    std::vector<float> l(kHop), r(kHop);
    for (std::uint32_t i = 0; i * kHop < kBlocks * 128; ++i) {
        OK(ms_source_set_pose(e, s, horizontal(i * 6.0, 2.0), {0, 0, 0}));
        OK(ms_engine_commit(e));
        ms_engine_render(e, l.data(), r.data(), kHop);
        left.insert(left.end(), l.begin(), l.end());
    }
    ms_engine_destroy(e);
    // Energy far above the tone, after a steep high pass, against total energy.
    std::vector<float> high = left;
    core::Biquad hp = core::designHighpass(4000.0f, 0.7071f, kRate);
    core::Biquad hp2 = hp, hp3 = hp, hp4 = hp;
    for (core::Biquad* q : {&hp, &hp2, &hp3, &hp4}) q->process(high.data(), static_cast<std::uint32_t>(high.size()));
    double highEnergy = 0.0, energy = 0.0, worstStep = 0.0;
    for (std::size_t i = 4800; i < left.size(); ++i) {
        highEnergy += static_cast<double>(high[i]) * high[i];
        energy += static_cast<double>(left[i]) * left[i];
        worstStep = std::max(worstStep, static_cast<double>(std::fabs(left[i] - left[i - 1])));
    }
    const double ratio = 10.0 * std::log10(highEnergy / energy);
    std::printf("  energy above 4 kHz with a 500 Hz tone orbiting: %.1f dB, largest sample step %.4f\n", ratio, worstStep);
    CHECK(ratio < -50.0, "%f dB", ratio);
}

void testVoiceOverflow() {
    gTest = "voice overflow";
    Scene scene;
    scene.build(gSphere, 4);
    for (int step = 0; step < 40; ++step) {
        scene.move(step);
        (void)renderFrames(scene.e, 1024);
    }
    ms_engine_stats stats;
    ms_engine_get_stats(scene.e, &stats);
    std::printf("  %u playing: %u convolved, %u through the shared mix\n", stats.sources_playing, stats.voices_direct, stats.voices_ambisonic);
    CHECK(stats.voices_direct <= 4, "%u", stats.voices_direct);
    CHECK(stats.voices_direct + stats.voices_ambisonic == stats.sources_playing, "voices lost");
    CHECK(stats.peak[0] > 0.0f, "silent");
    ms_engine_destroy(scene.e);
}

void testSchedulingAndFeatures() {
    gTest = "scheduling and features";
    std::vector<float> impulse(8, 0.0f);
    impulse[0] = 1.0f;
    // Sample accurate start: the output moves by exactly the requested frames.
    const auto startAt = [&](std::uint64_t frame) {
        ms_engine* e = makeEngine(gSphere);
        ms_buffer b;
        ms_source s;
        OK(ms_buffer_create(e, impulse.data(), impulse.size(), 1, kRate, &b));
        OK(ms_source_create(e, &s));
        const ms_source_params p = plainParams({0, 0, -1.5f});
        OK(ms_source_set_params(e, s, &p));
        OK(ms_source_set_buffer(e, s, b, 0, 0));
        OK(ms_source_play(e, s, frame));
        OK(ms_engine_commit(e));
        Render r = renderFrames(e, 4096);
        ms_engine_destroy(e);
        return r;
    };
    const Render at0 = startAt(0), at1000 = startAt(1000), at1037 = startAt(1037);
    const auto peakAt = [](const std::vector<float>& x) { return std::distance(x.begin(), std::max_element(x.begin(), x.end(), [](float a, float b) { return std::fabs(a) < std::fabs(b); })); };
    CHECK(peakAt(at1000.l) - peakAt(at0.l) == 1000, "moved by %ld", static_cast<long>(peakAt(at1000.l) - peakAt(at0.l)));
    CHECK(peakAt(at1037.l) - peakAt(at0.l) == 1037, "moved by %ld", static_cast<long>(peakAt(at1037.l) - peakAt(at0.l)));

    // The callback path with no Doppler matches the buffer path exactly.
    struct Feed {
        const std::vector<float>* data;
        std::size_t position;
        static uint32_t read(void* user, float* mono, uint32_t frames) {
            auto* f = static_cast<Feed*>(user);
            uint32_t i = 0;
            for (; i < frames && f->position < f->data->size(); ++i) mono[i] = (*f->data)[f->position++];
            return i;
        }
    } feed{&impulse, 0};
    ms_engine* e = makeEngine(gSphere);
    ms_source s;
    OK(ms_source_create(e, &s));
    const ms_source_params plain = plainParams({0, 0, -1.5f});
    OK(ms_source_set_params(e, s, &plain));
    OK(ms_source_set_callback(e, s, Feed::read, &feed));
    OK(ms_source_play(e, s, 0));
    OK(ms_engine_commit(e));
    const Render pulled = renderFrames(e, 4096);
    CHECK(mirrorError(at0.l, pulled.l) < -100.0, "callback differs %f dB", mirrorError(at0.l, pulled.l));
    ms_source_state state;
    OK(ms_source_get_state(e, s, &state));
    (void)renderFrames(e, 48000);
    OK(ms_source_get_state(e, s, &state));
    CHECK(state.playing == 0, "finished source still playing");
    ms_engine_destroy(e);

    // A bad callback cannot poison the output.
    struct Bad {
        static uint32_t read(void*, float* mono, uint32_t frames) {
            for (uint32_t i = 0; i < frames; ++i) mono[i] = i == 5 ? NAN : 0.1f;
            return frames;
        }
    };
    e = makeEngine(gSphere);
    OK(ms_source_create(e, &s));
    OK(ms_source_set_params(e, s, &plain));
    OK(ms_source_set_callback(e, s, Bad::read, nullptr));
    OK(ms_source_play(e, s, 0));
    OK(ms_engine_commit(e));
    const Render bad = renderFrames(e, 2048);
    bool finite = true;
    for (float v : bad.l) finite = finite && std::isfinite(v);
    ms_engine_stats stats;
    ms_engine_get_stats(e, &stats);
    CHECK(finite, "non-finite output");
    CHECK(stats.faults > 0, "fault not counted");
    ms_engine_destroy(e);

    // A reflection tap lands at its delay, from its side.
    e = makeEngine(gSphere, 128, 3);
    ms_buffer b;
    OK(ms_buffer_create(e, impulse.data(), impulse.size(), 1, kRate, &b));
    OK(ms_source_create(e, &s));
    ms_source_params silent = plainParams({0, 0, -1.5f});
    silent.gain = 1.0f;
    OK(ms_source_set_params(e, s, &silent));
    const ms_reflection_tap tap = {0.05f, {1.0f, 0.0f, 0.0f}, {1.0f, 1.0f, 1.0f}};
    OK(ms_source_set_reflections(e, s, &tap, 1));
    OK(ms_source_set_buffer(e, s, b, 0, 0));
    OK(ms_source_play(e, s, 0));
    OK(ms_engine_commit(e));
    const Render reflected = renderFrames(e, 9600);
    ms_engine_destroy(e);
    float tailL = 0.0f, tailR = 0.0f;
    std::size_t tailAt = 0;
    for (std::size_t i = 1200; i < reflected.r.size(); ++i) {
        if (std::fabs(reflected.r[i]) > tailR) {
            tailR = std::fabs(reflected.r[i]);
            tailAt = i;
        }
        tailL = std::max(tailL, std::fabs(reflected.l[i]));
    }
    const double tailDelay = static_cast<double>(tailAt) / kRate;
    std::printf("  reflection tap asked for 50.0 ms from the right: at %.2f ms, right %.3f left %.3f\n", tailDelay * 1000.0, tailR, tailL);
    CHECK(tailDelay > 0.049 && tailDelay < 0.0535, "%f s", tailDelay);
    CHECK(tailR > 1.5f * tailL, "not lateralised");

    // Bus gains multiply through the tree, and the limiter holds full scale.
    ms_engine_desc d;
    ms_engine_desc_default(&d);
    d.limiter = 1;
    d.ambisonic_order = 0;
    d.smoothing_seconds = 0.0f;
    OK(ms_engine_create(&d, &e));
    OK(ms_engine_set_hrtf(e, gSphere));
    std::vector<float> loud(9600, 40.0f);
    OK(ms_buffer_create(e, loud.data(), loud.size(), 1, kRate, &b));
    OK(ms_source_create(e, &s));
    OK(ms_source_set_params(e, s, &plain));
    OK(ms_source_set_buffer(e, s, b, 0, 0));
    OK(ms_source_play(e, s, 0));
    OK(ms_engine_commit(e));
    const Render limited = renderFrames(e, 4800);
    float peak = 0.0f;
    for (float v : limited.l) peak = std::max(peak, std::fabs(v));
    CHECK(peak <= 1.0f, "limiter let through %f", peak);
    ms_engine_destroy(e);
}

void testBlockSizesAgree() {
    gTest = "block sizes";
    std::vector<float> noise(4096);
    std::mt19937 rng(21);
    for (float& v : noise) v = std::uniform_real_distribution<float>(-0.5f, 0.5f)(rng);
    const auto run = [&](std::uint32_t block) {
        ms_engine* e = makeEngine(gSphere, block, 0);
        ms_buffer b;
        ms_source s;
        OK(ms_buffer_create(e, noise.data(), noise.size(), 1, kRate, &b));
        OK(ms_source_create(e, &s));
        const ms_source_params p = plainParams(horizontal(50.0, 1.5));
        OK(ms_source_set_params(e, s, &p));
        OK(ms_source_set_buffer(e, s, b, 0, 0));
        OK(ms_source_play(e, s, 0));
        OK(ms_engine_commit(e));
        Render r = renderFrames(e, 8192);
        ms_engine_destroy(e);
        return r;
    };
    const Render reference = run(128);
    for (std::uint32_t block : {64u, 256u, 512u}) {
        const double difference = mirrorError(reference.r, run(block).r);
        std::printf("  block %u against block 128: %.1f dB\n", block, difference);
        CHECK(difference < -100.0, "block %u differs by %f dB", block, difference);
    }
}

void testOtherRates() {
    gTest = "other rates";
    for (const float rate : {44100.0f, 96000.0f}) {
        ms_engine_desc d;
        ms_engine_desc_default(&d);
        d.sample_rate = rate;
        d.block_frames = rate > 50000.0f ? 256 : 128;
        d.ambisonic_order = 5;
        ms_engine* e = nullptr;
        OK(ms_engine_create(&d, &e));
        ms_hrtf* set = nullptr;
        OK(ms_hrtf_load(dataPath("sadie2_h3.mshrtf").c_str(), rate, &set));
        ms_hrtf_info info;
        ms_hrtf_get_info(set, &info);
        OK(ms_engine_set_hrtf(e, set));
        std::vector<float> stereo(4800);  // 48 kHz content, resampled on creation.
        for (std::size_t i = 0; i < stereo.size(); ++i) stereo[i] = 0.3f * std::sin(0.1f * static_cast<float>(i));
        ms_buffer b;
        ms_source s;
        OK(ms_buffer_create(e, stereo.data(), stereo.size(), 1, 48000.0f, &b));
        OK(ms_source_create(e, &s));
        const ms_source_params p = plainParams(horizontal(45.0, 1.2));
        OK(ms_source_set_params(e, s, &p));
        OK(ms_source_set_buffer(e, s, b, 0, 1));
        OK(ms_source_play(e, s, 0));
        OK(ms_engine_commit(e));
        std::vector<float> l(4096), r(4096);
        ms_engine_render(e, l.data(), r.data(), 4096);
        double energyL = 0.0, energyR = 0.0;
        for (std::size_t i = 1024; i < 4096; ++i) {
            energyL += static_cast<double>(l[i]) * l[i];
            energyR += static_cast<double>(r[i]) * r[i];
        }
        ms_engine_stats stats;
        ms_engine_get_stats(e, &stats);
        std::printf("  %.0f Hz: %u taps, right over left %.1f dB, %llu faults\n", static_cast<double>(rate), info.taps, 10.0 * std::log10(energyR / energyL), static_cast<unsigned long long>(stats.faults));
        CHECK(energyR > energyL && energyL > 0.0, "cue lost at %f Hz", static_cast<double>(rate));
        CHECK(stats.faults == 0, "faults");
        ms_hrtf_release(set);
        ms_engine_destroy(e);
    }
}

void testSourceShaping() {
    gTest = "source shaping";
    std::vector<float> impulse(8, 0.0f);
    impulse[0] = 1.0f;
    const auto render = [&](const std::function<void(ms_source_params&)>& shape, const std::function<void(ms_engine*)>& scene = {}) {
        ms_engine* e = makeEngine(gSphere, 128, 3);
        ms_buffer b;
        ms_source s;
        OK(ms_buffer_create(e, impulse.data(), impulse.size(), 1, kRate, &b));
        OK(ms_source_create(e, &s));
        ms_source_params p = plainParams({0, 0, -1.5f});
        shape(p);
        if (scene) scene(e);
        OK(ms_source_set_params(e, s, &p));
        OK(ms_source_set_buffer(e, s, b, 0, 0));
        OK(ms_source_play(e, s, 0));
        OK(ms_engine_commit(e));
        Render r = renderFrames(e, 4096);
        ms_engine_destroy(e);
        return r;
    };
    const Render plain = render([](ms_source_params&) {});
    // A cardioid facing away from the listener nulls; facing the listener passes.
    const Render facing = render([](ms_source_params& p) { p.directivity_alpha = 0.5f; p.forward = {0, 0, 1}; });
    const Render away = render([](ms_source_params& p) { p.directivity_alpha = 0.5f; p.forward = {0, 0, -1}; });
    CHECK(std::fabs(magnitudeDb(facing.l, 1000.0) - magnitudeDb(plain.l, 1000.0)) < 0.01, "cardioid on axis");
    CHECK(magnitudeDb(away.l, 1000.0) < magnitudeDb(plain.l, 1000.0) - 80.0, "cardioid rear %f", magnitudeDb(away.l, 1000.0));
    // Per-source filters do what their parameters say.
    const Render lowpassed = render([](ms_source_params& p) { p.lowpass_hz = 1000.0f; });
    CHECK(std::fabs(magnitudeDb(lowpassed.l, 1000.0) - magnitudeDb(plain.l, 1000.0) + 3.01) < 0.2, "lowpass corner");
    const Render highpassed = render([](ms_source_params& p) { p.highpass_hz = 1000.0f; });
    CHECK(std::fabs(magnitudeDb(highpassed.l, 1000.0) - magnitudeDb(plain.l, 1000.0) + 3.01) < 0.2, "highpass corner");
    const Render occluded = render([](ms_source_params& p) { p.eq_high_db = -12.0f; });
    CHECK(std::fabs(magnitudeDb(occluded.l, 12000.0) - magnitudeDb(plain.l, 12000.0) + 12.0) < 1.0, "high EQ");
    // Bus gains multiply through the tree.
    {
        ms_engine* e = makeEngine(gSphere, 128, 3);
        ms_bus a, b2;
        OK(ms_bus_create(e, MS_BUS_MASTER, &a));
        OK(ms_bus_create(e, a, &b2));
        OK(ms_bus_set_gain(e, a, 0.5f));
        OK(ms_bus_set_gain(e, b2, 0.5f));
        ms_buffer b;
        ms_source s;
        OK(ms_buffer_create(e, impulse.data(), impulse.size(), 1, kRate, &b));
        OK(ms_source_create(e, &s));
        ms_source_params p = plainParams({0, 0, -1.5f});
        p.bus = b2;
        OK(ms_source_set_params(e, s, &p));
        OK(ms_source_set_buffer(e, s, b, 0, 0));
        OK(ms_source_play(e, s, 0));
        OK(ms_engine_commit(e));
        const Render r = renderFrames(e, 4096);
        ms_engine_destroy(e);
        CHECK(std::fabs(magnitudeDb(r.l, 1000.0) - magnitudeDb(plain.l, 1000.0) + 12.04) < 0.05, "bus tree gain");
    }
    // Head tracking: a head turned 40 degrees right hears a source at 40 right as straight ahead.
    const float half = -20.0f * core::kPi / 180.0f;  // A negative rotation about +y turns a right-handed head to the right.
    const Render turned = render([](ms_source_params& p) { p.position = horizontal(40.0, 1.5); }, [&](ms_engine* e) { ms_head_set_orientation(e, {0, std::sin(half), 0, std::cos(half)}); });
    CHECK(mirrorError(plain.l, turned.l) < -40.0, "head turn differs by %f dB", mirrorError(plain.l, turned.l));
    // Pan mode ignores distance and the head.
    const Render panned = render([](ms_source_params& p) { p.mode = MS_SPATIAL_PAN; p.pan = 1.0f; p.position = {50, 0, 0}; });
    double left = 0.0, right = 0.0;
    for (std::size_t i = 0; i < panned.l.size(); ++i) {
        left += static_cast<double>(panned.l[i]) * panned.l[i];
        right += static_cast<double>(panned.r[i]) * panned.r[i];
    }
    CHECK(left < 1e-9 && std::fabs(right - 1.0) < 1e-3, "pan left %f right %f", left, right);
    // Swapping sets mid-stream stays clean, and a destroyed buffer's slot reuse is safe.
    ms_engine* e = makeEngine(gSphere, 128, 3);
    std::vector<float> tone(48000);
    for (std::size_t i = 0; i < tone.size(); ++i) tone[i] = 0.3f * std::sin(0.05f * static_cast<float>(i));
    ms_buffer b;
    ms_source s;
    OK(ms_buffer_create(e, tone.data(), tone.size(), 1, kRate, &b));
    OK(ms_source_create(e, &s));
    const ms_source_params plainSource = plainParams(horizontal(20.0, 1.5));
    OK(ms_source_set_params(e, s, &plainSource));
    OK(ms_source_set_buffer(e, s, b, 0, 1));
    OK(ms_source_play(e, s, 0));
    OK(ms_engine_commit(e));
    (void)renderFrames(e, 2048);
    OK(ms_engine_set_hrtf(e, gMeasured.empty() ? gSphere : gMeasured[0].second));
    OK(ms_engine_commit(e));
    (void)renderFrames(e, 2048);
    OK(ms_buffer_destroy(e, b));
    ms_buffer again;
    OK(ms_buffer_create(e, tone.data(), tone.size(), 1, kRate, &again));
    OK(ms_engine_commit(e));
    (void)renderFrames(e, 2048);
    OK(ms_engine_commit(e));
    ms_engine_stats stats;
    ms_engine_get_stats(e, &stats);
    CHECK(stats.faults == 0, "faults %llu", static_cast<unsigned long long>(stats.faults));
    ms_engine_destroy(e);
}

void testFilesAndSofa() {
    gTest = "files";
    const char* path = "/tmp/marathon_spatial_core_test.mshrtf";
    OK(ms_hrtf_save(gSphere, path));
    ms_hrtf* loaded = nullptr;
    OK(ms_hrtf_load(path, 0.0f, &loaded));
    ms_hrtf_info a, b;
    ms_hrtf_get_info(gSphere, &a);
    ms_hrtf_get_info(loaded, &b);
    CHECK(a.directions == b.directions && a.taps == b.taps && a.sample_rate == b.sample_rate, "info changed");
    const Render original = impulseThrough(gSphere, horizontal(33.0, 1.5)), reloaded = impulseThrough(loaded, horizontal(33.0, 1.5));
    CHECK(mirrorError(original.l, reloaded.l) < -120.0, "round trip changed the sound");
    ms_hrtf_release(loaded);

    // Resampling on load keeps the head.
    ms_hrtf* resampled = nullptr;
    OK(ms_hrtf_load(path, 44100.0f, &resampled));
    ms_hrtf_get_info(resampled, &b);
    CHECK(b.sample_rate == 44100.0f && b.directions == a.directions, "resampled info");
    ms_hrtf_release(resampled);

    // Damage is refused, not loaded.
    std::FILE* f = std::fopen(path, "r+b");
    std::fseek(f, 4096, SEEK_SET);
    std::fputc(0x5a, f);
    std::fclose(f);
    ms_hrtf* damaged = nullptr;
    CHECK(ms_hrtf_load(path, 0.0f, &damaged) == MS_ERR_BAD_FILE, "corrupt file accepted");
    CHECK(ms_hrtf_load("/nonexistent/x.mshrtf", 0.0f, &damaged) == MS_ERR_NOT_FOUND, "missing file");
    std::remove(path);

}

// A run of its own, so that ctest reports Skipped rather than Passed when there is no file: main returns 77 then.
void testSofaFile(const char* sofa) {
    gTest = "SOFA import";
    ms_hrtf_info b;
    {
        ms_hrtf* imported = nullptr;
        OK(ms_hrtf_import_sofa(sofa, nullptr, &imported));
        if (imported) {
            ms_hrtf_get_info(imported, &b);
            std::printf("  SOFA import of %s: %u directions at %.0f Hz, title '%s'\n", sofa, b.directions, static_cast<double>(b.sample_rate), b.name);
            const Render r = impulseThrough(imported, horizontal(90.0, 1.5));
            CHECK(onsetSeconds(r.l) > onsetSeconds(r.r), "SOFA ears crossed");
            ms_hrtf_release(imported);
        }
        // Truncated files must come back as errors.
        std::FILE* in = std::fopen(sofa, "rb");
        std::fseek(in, 0, SEEK_END);
        const long size = std::ftell(in);
        std::fseek(in, 0, SEEK_SET);
        std::vector<char> bytes(static_cast<std::size_t>(size));
        (void)std::fread(bytes.data(), 1, bytes.size(), in);
        std::fclose(in);
        for (const double keep : {0.001, 0.1, 0.5, 0.9}) {
            const char* cut = "/tmp/marathon_spatial_core_test.sofa";
            std::FILE* out = std::fopen(cut, "wb");
            std::fwrite(bytes.data(), 1, static_cast<std::size_t>(static_cast<double>(size) * keep), out);
            std::fclose(out);
            ms_hrtf* partial = nullptr;
            CHECK(ms_hrtf_import_sofa(cut, nullptr, &partial) != MS_OK, "truncated SOFA accepted at %f", keep);
            std::remove(cut);
        }
    }
}

// A source whose direct sound the distance curve silences, while its dry signal, and so its reflections, keep playing.
ms_source_params tapOnlyParams(ms_vec3 position = {0, 0, -10}) {
    ms_source_params p = plainParams(position);
    p.distance_model = MS_DISTANCE_LINEAR;
    p.ref_distance = 1.0f;
    p.max_distance = 2.0f;
    p.rolloff = 1.0f;
    return p;
}

std::vector<float> noiseOf(std::size_t frames, std::uint32_t seed = 1) {
    std::mt19937 random(seed);
    std::uniform_real_distribution<float> sample(-0.5f, 0.5f);
    std::vector<float> x(frames);
    for (float& v : x) v = sample(random);
    return x;
}

std::vector<float> sineOf(std::size_t frames, double hz, double rate = kRate) {
    std::vector<float> x(frames);
    for (std::size_t i = 0; i < frames; ++i) x[i] = 0.3f * static_cast<float>(std::sin(2.0 * M_PI * hz * static_cast<double>(i) / rate));
    return x;
}

double rmsDb(const std::vector<float>& a, const std::vector<float>& b, std::size_t from) {
    double sum = 0.0;
    for (std::size_t i = from; i < a.size(); ++i) sum += static_cast<double>(a[i]) * a[i] + static_cast<double>(b[i]) * b[i];
    return 10.0 * std::log10(sum / static_cast<double>(2 * (a.size() - from)) + 1e-30);
}

double largestStep(const std::vector<float>& x, std::size_t from, std::size_t to) {
    double step = 0.0;
    for (std::size_t i = std::max<std::size_t>(from, 1); i < to; ++i) step = std::max(step, std::fabs(static_cast<double>(x[i]) - x[i - 1]));
    return step;
}

// Finding: ambisonic content went through the master's gain twice, once in the bus product and once on the master bus.
void testMasterGain() {
    gTest = "master gain";
    const auto level = [](std::uint32_t voices, float master, float child) {
        ms_engine* e = makeEngine(gSphere, 128, 3, 0.0f, voices);
        const std::vector<float> noise = noiseOf(48000);
        ms_buffer b;
        ms_source s;
        ms_bus bus;
        OK(ms_buffer_create(e, noise.data(), noise.size(), 1, kRate, &b));
        OK(ms_source_create(e, &s));
        OK(ms_bus_create(e, MS_BUS_MASTER, &bus));
        OK(ms_bus_set_gain(e, MS_BUS_MASTER, master));
        OK(ms_bus_set_gain(e, bus, child));
        ms_source_params p = plainParams(horizontal(30.0, 2.0));
        p.bus = bus;
        OK(ms_source_set_params(e, s, &p));
        OK(ms_source_set_buffer(e, s, b, 0, 1));
        OK(ms_source_play(e, s, 0));
        OK(ms_engine_commit(e));
        const Render r = renderFrames(e, 16384);
        ms_engine_destroy(e);
        return rmsDb(r.l, r.r, 4096);
    };
    const double overflow = level(0, 1.0f, 1.0f), direct = level(8, 1.0f, 1.0f);
    const double overflowMaster = level(0, 0.25f, 1.0f) - overflow, directMaster = level(8, 0.25f, 1.0f) - direct;
    const double overflowBoth = level(0, 0.5f, 0.5f) - overflow;
    std::printf("  master -12 dB: overflow %.2f dB, direct %.2f dB; master and bus -6 dB each: overflow %.2f dB\n", overflowMaster, directMaster, overflowBoth);
    CHECK(std::fabs(overflowMaster + 12.04) < 0.1, "overflow voices at %.2f dB", overflowMaster);
    CHECK(std::fabs(directMaster + 12.04) < 0.1, "direct voices at %.2f dB", directMaster);
    CHECK(std::fabs(overflowBoth + 12.04) < 0.1, "bus and master at %.2f dB", overflowBoth);
}

// Finding: reflections and beds reach the master without their bus, so ducking or muting a group left them audible.
void testBusesReachEverything() {
    gTest = "buses reach everything";
    const auto tapLevel = [](float gain) {
        ms_engine* e = makeEngine(gSphere, 128, 3);
        const std::vector<float> noise = noiseOf(48000);
        ms_buffer b;
        ms_source s;
        ms_bus group, child;
        OK(ms_buffer_create(e, noise.data(), noise.size(), 1, kRate, &b));
        OK(ms_source_create(e, &s));
        OK(ms_bus_create(e, MS_BUS_MASTER, &group));
        OK(ms_bus_create(e, group, &child));
        OK(ms_bus_set_gain(e, group, gain));
        ms_source_params p = tapOnlyParams();
        p.bus = child;
        OK(ms_source_set_params(e, s, &p));
        const ms_reflection_tap tap = {0.01f, {1.0f, 0.0f, 0.0f}, {1.0f, 1.0f, 1.0f}};
        OK(ms_source_set_reflections(e, s, &tap, 1));
        OK(ms_source_set_buffer(e, s, b, 0, 1));
        OK(ms_source_play(e, s, 0));
        OK(ms_engine_commit(e));
        const Render r = renderFrames(e, 16384);
        ms_engine_destroy(e);
        return rmsDb(r.l, r.r, 4096);
    };
    const double tapFull = tapLevel(1.0f), tapHalf = tapLevel(0.5f), tapMuted = tapLevel(0.0f);
    std::printf("  a tap through a group: full %.1f dB, group -6 dB %.2f dB, muted %.1f dB\n", tapFull, tapHalf - tapFull, tapMuted);
    CHECK(tapFull > -80.0, "the tap made no sound");
    CHECK(std::fabs(tapHalf - tapFull + 6.02) < 0.1, "%.2f dB", tapHalf - tapFull);
    CHECK(tapMuted < -200.0, "a muted group still sounds at %.1f dB", tapMuted);

    const auto bedLevel = [](float gain) {
        ms_engine* e = makeEngine(gSphere, 128, 1);
        const std::vector<float> noise = noiseOf(48000);
        std::vector<float> channels(noise.size() * 4, 0.0f);
        for (std::size_t i = 0; i < noise.size(); ++i) channels[i * 4] = noise[i];
        ms_buffer b;
        ms_bed bed;
        ms_bus group;
        OK(ms_buffer_create(e, channels.data(), noise.size(), 4, kRate, &b));
        OK(ms_bed_create(e, &bed));
        OK(ms_bus_create(e, MS_BUS_MASTER, &group));
        OK(ms_bus_set_gain(e, group, gain));
        OK(ms_bed_set_bus(e, bed, group));
        OK(ms_bed_set_buffer(e, bed, b, 1));
        OK(ms_bed_play(e, bed, 0));
        OK(ms_engine_commit(e));
        const Render r = renderFrames(e, 16384);
        ms_engine_destroy(e);
        return rmsDb(r.l, r.r, 4096);
    };
    const double bedFull = bedLevel(1.0f), bedHalf = bedLevel(0.5f);
    std::printf("  a bed on a bus at -6 dB: %.2f dB\n", bedHalf - bedFull);
    CHECK(std::fabs(bedHalf - bedFull + 6.02) < 0.1, "%.2f dB", bedHalf - bedFull);
    ms_engine* e = makeEngine(gSphere, 128, 1);
    ms_bed bed;
    OK(ms_bed_create(e, &bed));
    CHECK(ms_bed_set_bus(e, bed, 0x7fff0002u) == MS_ERR_NOT_FOUND, "a dead bus was accepted");
    ms_engine_destroy(e);
}

struct FirstPull {
    std::atomic<std::uint64_t>* clock;
    std::atomic<std::uint64_t> block{~0ull};
};
std::uint32_t firstPullCallback(void* user, float* mono, std::uint32_t frames) {
    auto* f = static_cast<FirstPull*>(user);
    std::uint64_t none = ~0ull;
    f->block.compare_exchange_strong(none, f->clock->load(std::memory_order_relaxed), std::memory_order_relaxed);
    std::memset(mono, 0, frames * sizeof(float));
    return frames;
}
void countBlocks(void* user, std::uint32_t, std::uint64_t) { static_cast<std::atomic<std::uint64_t>*>(user)->fetch_add(1, std::memory_order_relaxed); }

// Finding: commands were published one at a time, so a stereo pair could lock a block apart for its whole life.
void testCommitIsOneBatch() {
    gTest = "commit is one batch";
    core::SpscRing<std::uint32_t> ring;
    ring.init(4);
    const std::uint32_t three[3] = {1, 2, 3};
    CHECK(ring.pushAll(three, 3), "three fit in four");
    CHECK(!ring.pushAll(three, 2) && ring.pending() == 3, "a batch that does not fit must leave nothing behind");
    CHECK(ring.pushAll(three, 1) && ring.pending() == 4, "the last slot is usable");

    // One slot short of room: the pair must stay together, whatever the ring does with it.
    std::atomic<std::uint64_t> clock{0};
    ms_engine* e = makeEngine(gSphere, 128, 3);
    FirstPull left{&clock}, right{&clock};
    ms_source l, r, filler;
    OK(ms_source_create(e, &l));
    OK(ms_source_create(e, &r));
    OK(ms_source_create(e, &filler));
    OK(ms_engine_set_block_callback(e, countBlocks, &clock));
    OK(ms_source_set_callback(e, l, firstPullCallback, &left));
    OK(ms_source_set_callback(e, r, firstPullCallback, &right));
    OK(ms_engine_commit(e));
    (void)renderFrames(e, 128);
    const std::uint32_t room = 8192;
    for (std::uint32_t i = 0; i + 1 < room; ++i) OK(ms_source_seek(e, filler, 0.0));
    OK(ms_engine_commit(e));
    ms_engine_stats stats;
    ms_engine_get_stats(e, &stats);
    CHECK(stats.commands_pending == room - 1, "%u pending", stats.commands_pending);
    OK(ms_source_play(e, l, 0));
    OK(ms_source_play(e, r, 0));
    OK(ms_engine_commit(e));
    for (int i = 0; i < 4; ++i) {
        (void)renderFrames(e, 128);
        OK(ms_engine_commit(e));
    }
    std::printf("  with one slot free the pair started in blocks %llu and %llu\n", static_cast<unsigned long long>(left.block.load()), static_cast<unsigned long long>(right.block.load()));
    CHECK(left.block.load() != ~0ull && left.block.load() == right.block.load(), "split: %llu and %llu", static_cast<unsigned long long>(left.block.load()), static_cast<unsigned long long>(right.block.load()));
    ms_engine_destroy(e);

    // The same with a render thread racing the commits.
    std::atomic<std::uint64_t> raceClock{0};
    e = makeEngine(gSphere, 128, 3);
    FirstPull a{&raceClock}, b{&raceClock};
    OK(ms_source_create(e, &l));
    OK(ms_source_create(e, &r));
    OK(ms_engine_set_block_callback(e, countBlocks, &raceClock));
    OK(ms_source_set_callback(e, l, firstPullCallback, &a));
    OK(ms_source_set_callback(e, r, firstPullCallback, &b));
    OK(ms_engine_commit(e));
    std::atomic<bool> running{true};
    std::thread render([&] {
        std::vector<float> x(128), y(128);
        while (running.load(std::memory_order_relaxed)) ms_engine_render(e, x.data(), y.data(), 128);
    });
    int split = 0, rounds = 0;
    for (; rounds < 300; ++rounds) {
        a.block.store(~0ull, std::memory_order_relaxed);
        b.block.store(~0ull, std::memory_order_relaxed);
        OK(ms_source_play(e, l, 0));
        OK(ms_source_play(e, r, 0));
        OK(ms_engine_commit(e));
        while (a.block.load(std::memory_order_relaxed) == ~0ull || b.block.load(std::memory_order_relaxed) == ~0ull) std::this_thread::yield();
        split += a.block.load(std::memory_order_relaxed) != b.block.load(std::memory_order_relaxed);
        OK(ms_source_stop(e, l));
        OK(ms_source_stop(e, r));
        OK(ms_engine_commit(e));
        // Both idle, or a late pull of this round would be taken for the next round's start.
        ms_source_state stateL, stateR;
        do {
            std::this_thread::yield();
            OK(ms_engine_commit(e));
            OK(ms_source_get_state(e, l, &stateL));
            OK(ms_source_get_state(e, r, &stateR));
        } while (stateL.playing || stateR.playing);
    }
    running = false;
    render.join();
    std::printf("  %d pairs started while rendering: %d split\n", rounds, split);
    CHECK(split == 0, "%d pairs split", split);
    ms_engine_destroy(e);
}

struct BusEnergy {
    double sum = 0.0;
};
void measureBus(void* user, float* left, float*, std::uint32_t frames, std::uint64_t) {
    auto* m = static_cast<BusEnergy*>(user);
    for (std::uint32_t i = 0; i < frames; ++i) m->sum += static_cast<double>(left[i]) * left[i];
}

// Finding: a send was taken after the listener path's filters, so a source behind a door barely excited its own room.
void testSendsHearTheSource() {
    gTest = "sends hear the source";
    const auto sent = [](bool occluded) {
        ms_engine* e = makeEngine(gSphere, 128, 3);
        const std::vector<float> noise = noiseOf(48000);
        ms_buffer b;
        ms_source s;
        ms_bus room;
        BusEnergy energy;
        OK(ms_buffer_create(e, noise.data(), noise.size(), 1, kRate, &b));
        OK(ms_source_create(e, &s));
        OK(ms_bus_create(e, MS_BUS_MASTER, &room));
        OK(ms_bus_set_processor(e, room, measureBus, &energy));
        ms_source_params p = plainParams(horizontal(30.0, 2.0));
        p.gain = 0.5f;
        p.send_bus[0] = room;
        p.send_gain[0] = 0.5f;
        if (occluded) {
            p.lowpass_hz = 300.0f;
            p.eq_low_db = -24.0f;
            p.eq_mid_db = -24.0f;
            p.eq_high_db = -24.0f;
            p.air_absorption = 1.0f;
        }
        OK(ms_source_set_params(e, s, &p));
        OK(ms_source_set_buffer(e, s, b, 0, 1));
        OK(ms_source_play(e, s, 0));
        OK(ms_engine_commit(e));
        (void)renderFrames(e, 16384);
        ms_engine_destroy(e);
        return 10.0 * std::log10(energy.sum);
    };
    const double clear = sent(false), behind = sent(true);
    std::printf("  send energy behind a door relative to a clear line: %.3f dB\n", behind - clear);
    CHECK(std::fabs(behind - clear) < 0.01, "%.2f dB", behind - clear);
}

void poisonBus(void*, float* left, float* right, std::uint32_t, std::uint64_t) {
    left[3] = NAN;
    right[5] = INFINITY;
}

// Finding: a bad value reached a shared bus, and the guard at the end then silenced the whole mix.
void testBadValuesStayLocal() {
    gTest = "bad values stay local";
    ms_engine* e = makeEngine(gSphere, 128, 3);
    ms_source s;
    OK(ms_source_create(e, &s));
    float ms_source_params::*const fields[] = {&ms_source_params::gain, &ms_source_params::pitch, &ms_source_params::pan, &ms_source_params::ref_distance,
                                               &ms_source_params::max_distance, &ms_source_params::rolloff, &ms_source_params::doppler,
                                               &ms_source_params::air_absorption, &ms_source_params::near_field, &ms_source_params::directivity_alpha,
                                               &ms_source_params::directivity_sharpness, &ms_source_params::directivity_hf_db, &ms_source_params::highpass_hz,
                                               &ms_source_params::lowpass_hz, &ms_source_params::eq_low_db, &ms_source_params::eq_mid_db,
                                               &ms_source_params::eq_high_db};
    int accepted = 0;
    for (auto field : fields)
        for (const float bad : {NAN, INFINITY}) {
            ms_source_params p = plainParams({0, 0, -1});
            p.*field = bad;
            accepted += ms_source_set_params(e, s, &p) == MS_OK;
        }
    ms_source_params sendBad = plainParams({0, 0, -1});
    sendBad.send_gain[1] = NAN;
    accepted += ms_source_set_params(e, s, &sendBad) == MS_OK;
    ms_source_params negative = plainParams({0, 0, -1});
    negative.lowpass_hz = -1.0f;
    accepted += ms_source_set_params(e, s, &negative) == MS_OK;
    for (int i = 0; i < 6; ++i) {
        ms_reflection_tap tap = {0.01f, {1, 0, 0}, {1, 1, 1}};
        if (i == 0) tap.delay_seconds = NAN;
        if (i == 1) tap.delay_seconds = -0.01f;
        if (i == 2) tap.direction.y = INFINITY;
        if (i >= 3) tap.gain[i - 3] = NAN;
        accepted += ms_source_set_reflections(e, s, &tap, 1) == MS_OK;
    }
    ms_bed bed;
    OK(ms_bed_create(e, &bed));
    accepted += ms_bed_set_orientation(e, bed, ms_quat{NAN, 0, 0, 1}, 0) == MS_OK;
    CHECK(accepted == 0, "%d bad values accepted", accepted);
    ms_engine_destroy(e);

    // One source overflows, the other plays on. Every block of the mix must keep sounding.
    const auto clean = [](ms_engine* engine, ms_source& source, ms_buffer& buffer) {
        const std::vector<float> noise = noiseOf(48000, 3);
        OK(ms_buffer_create(engine, noise.data(), noise.size(), 1, kRate, &buffer));
        OK(ms_source_create(engine, &source));
        const ms_source_params p = plainParams(horizontal(-30.0, 2.0));
        OK(ms_source_set_params(engine, source, &p));
        OK(ms_source_set_buffer(engine, source, buffer, 0, 1));
        OK(ms_source_play(engine, source, 0));
    };
    const auto silentBlocks = [](ms_engine* engine) {
        int silent = 0;
        (void)renderFrames(engine, 2048);
        for (int block = 0; block < 64; ++block) {
            const Render r = renderFrames(engine, 128);
            double sum = 0.0;
            for (std::size_t i = 0; i < 128; ++i) sum += static_cast<double>(r.l[i]) * r.l[i];
            silent += !(sum > 0.0);
        }
        return silent;
    };
    e = makeEngine(gSphere, 128, 3);
    ms_source good, huge;
    ms_buffer content, loud;
    clean(e, good, content);
    const std::vector<float> noise = noiseOf(48000, 4);
    OK(ms_buffer_create(e, noise.data(), noise.size(), 1, kRate, &loud));
    OK(ms_source_create(e, &huge));
    ms_source_params overflowing = plainParams(horizontal(30.0, 2.0));
    overflowing.gain = 1e30f;
    OK(ms_source_set_params(e, huge, &overflowing));
    OK(ms_source_set_buffer(e, huge, loud, 0, 1));
    OK(ms_source_play(e, huge, 0));
    OK(ms_engine_commit(e));
    int silent = silentBlocks(e);
    ms_engine_stats stats;
    ms_engine_get_stats(e, &stats);
    std::printf("  overflowing source: %d of 64 mix blocks silent, %llu faults\n", silent, static_cast<unsigned long long>(stats.faults));
    CHECK(silent == 0, "the mix went silent for %d blocks", silent);
    CHECK(stats.faults > 0, "fault not counted");
    ms_engine_destroy(e);

    // Finite taps whose sum overflows lose their own reflections, not the shared ambisonic mix.
    e = makeEngine(gSphere, 128, 3);
    clean(e, good, content);
    ms_reflection_tap opposing[2] = {};
    for (int i = 0; i < 2; ++i) {
        opposing[i].delay_seconds = 0.01f;
        opposing[i].direction = {i ? -1.0f : 1.0f, 0.0f, 0.0f};
        opposing[i].gain[0] = opposing[i].gain[1] = opposing[i].gain[2] = 3e38f;
    }
    OK(ms_source_set_reflections(e, good, opposing, 2));
    OK(ms_engine_commit(e));
    silent = silentBlocks(e);
    ms_engine_get_stats(e, &stats);
    std::printf("  overflowing taps: %d of 64 mix blocks silent, %llu faults\n", silent, static_cast<unsigned long long>(stats.faults));
    CHECK(silent == 0, "the mix went silent for %d blocks", silent);
    CHECK(stats.faults > 0, "fault not counted");
    ms_engine_destroy(e);

    // A processor that writes NaN loses its own bus, not the master.
    e = makeEngine(gSphere, 128, 3);
    clean(e, good, content);
    ms_bus broken;
    OK(ms_bus_create(e, MS_BUS_MASTER, &broken));
    OK(ms_bus_set_processor(e, broken, poisonBus, nullptr));
    OK(ms_engine_commit(e));
    silent = silentBlocks(e);
    std::printf("  poisoned bus: %d of 64 mix blocks silent\n", silent);
    CHECK(silent == 0, "the mix went silent for %d blocks", silent);
    ms_engine_destroy(e);

    // So does a bed with a bad sample, which would otherwise sit in the decoder's history.
    e = makeEngine(gSphere, 128, 3);
    clean(e, good, content);
    std::vector<float> poisoned(16 * 4800, 0.1f);
    poisoned[16 * 2600] = NAN;
    OK(ms_buffer_create(e, poisoned.data(), 4800, 16, kRate, &loud));
    OK(ms_bed_create(e, &bed));
    OK(ms_bed_set_buffer(e, bed, loud, 1));
    OK(ms_bed_play(e, bed, 0));
    OK(ms_engine_commit(e));
    silent = silentBlocks(e);
    std::printf("  poisoned bed: %d of 64 mix blocks silent\n", silent);
    CHECK(silent == 0, "the mix went silent for %d blocks", silent);
    ms_engine_destroy(e);

    // A finite EQ can still design a filter with non-finite coefficients. That source goes quiet, and nothing else does.
    e = makeEngine(gSphere, 128, 3);
    clean(e, good, content);
    OK(ms_source_create(e, &huge));
    ms_source_params extreme = plainParams(horizontal(30.0, 2.0));
    extreme.eq_low_db = 1e6f;
    OK(ms_source_set_params(e, huge, &extreme));
    OK(ms_source_set_buffer(e, huge, content, 0, 1));
    OK(ms_source_play(e, huge, 0));
    OK(ms_engine_commit(e));
    silent = silentBlocks(e);
    ms_engine_get_stats(e, &stats);
    std::printf("  extreme EQ: %d of 64 mix blocks silent, %llu faults\n", silent, static_cast<unsigned long long>(stats.faults));
    CHECK(silent == 0, "the mix went silent for %d blocks", silent);
    CHECK(stats.faults > 0, "this proves nothing unless the filter really went bad");
    ms_engine_destroy(e);
}

// Finding: losing a direct voice threw away its filter tail and every reflection still in its delay lines, and two sources
// of about the same level traded the voice at every ranking.
void testDemotionKeepsTails() {
    gTest = "demotion keeps tails";
    const auto run = [](bool demote, bool tapOnly, std::uint32_t impulseAt, std::uint32_t frames, std::uint32_t* directAfter) {
        ms_engine* e = makeEngine(gSphere, 128, 3, 0.0f, 1);
        std::vector<float> content(frames + 2048, 0.0f), silence(4800, 0.0f);
        content[impulseAt] = 1.0f;
        ms_buffer b, quiet;
        ms_source a, other;
        OK(ms_buffer_create(e, content.data(), content.size(), 1, kRate, &b));
        OK(ms_buffer_create(e, silence.data(), silence.size(), 1, kRate, &quiet));
        OK(ms_source_create(e, &a));
        OK(ms_source_create(e, &other));
        const ms_source_params p = tapOnly ? tapOnlyParams() : plainParams(horizontal(30.0, 2.0));
        OK(ms_source_set_params(e, a, &p));
        const ms_reflection_tap tap = {0.06f, {1, 0, 0}, {1, 1, 1}};
        if (tapOnly) OK(ms_source_set_reflections(e, a, &tap, 1));
        OK(ms_source_set_buffer(e, a, b, 0, 0));
        OK(ms_source_play(e, a, 0));
        OK(ms_engine_commit(e));
        Render out = renderFrames(e, 512);
        if (demote) {  // Silent, and far louder by its gain, so it only takes the voice away.
            ms_source_params loud = plainParams(horizontal(-30.0, 2.0));
            loud.gain = 100.0f;
            OK(ms_source_set_params(e, other, &loud));
            OK(ms_source_set_buffer(e, other, quiet, 0, 1));
            OK(ms_source_play(e, other, 0));
            OK(ms_engine_commit(e));
        }
        const Render rest = renderFrames(e, frames - 512);
        out.l.insert(out.l.end(), rest.l.begin(), rest.l.end());
        out.r.insert(out.r.end(), rest.r.begin(), rest.r.end());
        ms_source_state state;
        OK(ms_source_get_state(e, a, &state));
        if (directAfter) *directAfter = state.direct;
        ms_engine_destroy(e);
        return out;
    };
    // The reflection of a sample heard before the demotion, 60 ms later.
    std::uint32_t directAfter = 1;
    const Render kept = run(false, true, 256, 9600, nullptr), demoted = run(true, true, 256, 9600, &directAfter);
    const auto energy = [](const Render& r) {
        double sum = 0.0;
        for (std::size_t i = 2048; i < r.l.size(); ++i) sum += static_cast<double>(r.l[i]) * r.l[i] + static_cast<double>(r.r[i]) * r.r[i];
        return sum;
    };
    const double ratio = 10.0 * std::log10(energy(demoted) / energy(kept) + 1e-30);
    std::printf("  a reflection in flight when its voice is demoted: %.3f dB against no demotion\n", ratio);
    CHECK(directAfter == 0, "the voice was never demoted, so this proves nothing");
    CHECK(std::fabs(ratio) < 0.05, "%.2f dB", ratio);

    // The filter tail: nothing goes in once the voice is demoted, so what rings out must be identical.
    const Render tailKept = run(false, false, 1023, 2048, nullptr), tailDemoted = run(true, false, 1023, 2048, &directAfter);
    double worst = 0.0;
    for (std::size_t i = 0; i < tailKept.l.size(); ++i) {
        worst = std::max(worst, std::fabs(static_cast<double>(tailKept.l[i]) - tailDemoted.l[i]));
        worst = std::max(worst, std::fabs(static_cast<double>(tailKept.r[i]) - tailDemoted.r[i]));
    }
    std::printf("  filter tail across a demotion: worst difference %.3g\n", worst);
    CHECK(directAfter == 0, "the voice was never demoted, so this proves nothing");
    CHECK(worst < 1e-5, "tail cut: %g", worst);

    // A held voice is kept against something only a little louder, and given up to something clearly louder.
    ms_engine* e = makeEngine(gSphere, 128, 3, 0.0f, 1);
    const std::vector<float> noise = noiseOf(48000);
    ms_buffer b;
    ms_source a, other;
    OK(ms_buffer_create(e, noise.data(), noise.size(), 1, kRate, &b));
    OK(ms_source_create(e, &a));
    OK(ms_source_create(e, &other));
    const ms_source_params first = plainParams(horizontal(30.0, 2.0));
    ms_source_params second = plainParams(horizontal(-30.0, 2.0));
    second.gain = 1.25f;
    OK(ms_source_set_params(e, a, &first));
    OK(ms_source_set_buffer(e, a, b, 0, 1));
    OK(ms_source_play(e, a, 0));
    OK(ms_engine_commit(e));
    (void)renderFrames(e, 4096);
    OK(ms_source_set_params(e, other, &second));
    OK(ms_source_set_buffer(e, other, b, 0, 1));
    OK(ms_source_play(e, other, 0));
    OK(ms_engine_commit(e));
    (void)renderFrames(e, 8192);
    ms_source_state stateA, stateB;
    OK(ms_source_get_state(e, a, &stateA));
    OK(ms_source_get_state(e, other, &stateB));
    CHECK(stateA.direct == 1 && stateB.direct == 0, "traded at +2 dB: %u %u", stateA.direct, stateB.direct);
    second.gain = 2.0f;
    OK(ms_source_set_params(e, other, &second));
    OK(ms_engine_commit(e));
    (void)renderFrames(e, 16384);
    OK(ms_source_get_state(e, a, &stateA));
    OK(ms_source_get_state(e, other, &stateB));
    CHECK(stateA.direct == 0 && stateB.direct == 1, "kept at +6 dB: %u %u", stateA.direct, stateB.direct);
    ms_engine_destroy(e);
}

// Finding: the reflections' copy of the signal jumped on a gain change and on a stop, which is a click in every tap.
void testDrySignalRamps() {
    gTest = "dry signal ramps";
    constexpr std::uint32_t kChange = 8192;
    enum class Change { None, ContentStep, Gain, Stop };
    const auto run = [](Change change) {
        ms_engine* e = makeEngine(gSphere, 128, 3);
        std::vector<float> content(48000, 0.25f);
        if (change == Change::ContentStep) std::fill(content.begin() + kChange, content.end(), 0.0f);
        ms_buffer b;
        ms_source s;
        OK(ms_buffer_create(e, content.data(), content.size(), 1, kRate, &b));
        OK(ms_source_create(e, &s));
        ms_source_params p = tapOnlyParams();
        OK(ms_source_set_params(e, s, &p));
        const ms_reflection_tap tap = {0.01f, {1, 0, 0}, {1, 1, 1}};
        OK(ms_source_set_reflections(e, s, &tap, 1));
        OK(ms_source_set_buffer(e, s, b, 0, 1));
        OK(ms_source_play(e, s, 0));
        OK(ms_engine_commit(e));
        Render out = renderFrames(e, kChange);
        if (change == Change::Gain) {
            p.gain = 0.0f;
            OK(ms_source_set_params(e, s, &p));
        } else if (change == Change::Stop) {
            OK(ms_source_stop(e, s));
        }
        OK(ms_engine_commit(e));
        const Render rest = renderFrames(e, 4096);
        out.r.insert(out.r.end(), rest.r.begin(), rest.r.end());
        ms_engine_destroy(e);
        return largestStep(out.r, kChange, out.r.size());
    };
    const double abrupt = run(Change::ContentStep), gain = run(Change::Gain), stopped = run(Change::Stop);
    std::printf("  largest step in the tap: content cut %.4f, gain to zero %.4f, stop %.4f\n", abrupt, gain, stopped);
    CHECK(abrupt > 0.0, "the reference made no step");
    CHECK(gain < 0.5 * abrupt, "gain change clicks: %f against %f", gain, abrupt);
    CHECK(stopped < 0.5 * abrupt, "stop clicks: %f against %f", stopped, abrupt);
}

std::uint32_t countFrames(void* user, float* mono, std::uint32_t frames) {
    *static_cast<std::uint64_t*>(user) += frames;
    std::memset(mono, 0, frames * sizeof(float));
    return frames;
}

// Finding: pitch stopped at 0.25, so 16 kHz speech on a 96 kHz device ran 1.5 times fast.
void testLowPitch() {
    gTest = "low pitch";
    ms_engine* e = makeEngine(gSphere, 128, 3);
    const std::vector<float> noise = noiseOf(48000);
    ms_buffer b;
    ms_source buffered, pulled;
    OK(ms_buffer_create(e, noise.data(), noise.size(), 1, kRate, &b));
    OK(ms_source_create(e, &buffered));
    OK(ms_source_create(e, &pulled));
    ms_source_params p = plainParams({0, 0, -1});
    p.pitch = 1.0f / 6.0f;
    OK(ms_source_set_params(e, buffered, &p));
    OK(ms_source_set_params(e, pulled, &p));
    std::uint64_t frames = 0;
    OK(ms_source_set_buffer(e, buffered, b, 0, 1));
    OK(ms_source_set_callback(e, pulled, countFrames, &frames));
    OK(ms_source_play(e, buffered, 0));
    OK(ms_source_play(e, pulled, 0));
    OK(ms_engine_commit(e));
    (void)renderFrames(e, 9600);
    ms_source_state state;
    OK(ms_source_get_state(e, buffered, &state));
    std::printf("  9600 frames at pitch 1/6: buffer advanced %.0f, callback gave %llu\n", state.position_frames, static_cast<unsigned long long>(frames));
    CHECK(std::fabs(state.position_frames - 1600.0) < 2.0, "%f", state.position_frames);
    CHECK(frames >= 1600 && frames < 1640, "%llu", static_cast<unsigned long long>(frames));
    ms_engine_destroy(e);
}

// Finding: on EALREADY the device thread asked again on every cycle.
void testWorkgroupJoinOnce() {
    gTest = "workgroup join";
    using Result = core::WorkgroupJoin::Result;
    for (const int answer : {EALREADY, EINVAL}) {
        core::WorkgroupJoin join;
        int asked = 0;
        Result first = Result::NotAttempted;
        for (int cycle = 0; cycle < 100; ++cycle) {
            const Result r = join.step([&] {
                ++asked;
                return answer;
            });
            if (cycle == 0) first = r;
        }
        CHECK(asked == 1, "asked %d times", asked);
        CHECK(first == (answer == EALREADY ? Result::AlreadyMember : Result::Failed), "wrong answer");
        join.reset();
        (void)join.step([&] {
            ++asked;
            return 0;
        });
        CHECK(asked == 2, "a restarted device thread must ask again");
    }
    core::WorkgroupJoin join;
    CHECK(join.step([] { return 0; }) == Result::Joined, "a good join");
}

// Six directions of one impulse each. What the build does with them is then known exactly.
hrtf::RawHrirSet rawImpulses(float rate, std::uint32_t taps) {
    hrtf::RawHrirSet raw;
    raw.taps = taps;
    raw.sampleRate = rate;
    raw.distance = 1.5f;
    raw.directions = {{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}};
    raw.samples.assign(raw.directions.size() * 2 * taps, 0.0f);
    for (std::size_t d = 0; d < raw.directions.size() * 2; ++d) raw.samples[d * taps + 32] = 1.0f;
    return raw;
}

// Finding of this round: Data.Delay was dropped on import, and a file that keeps its ITD there came in with none.
void testStoredDelaysSurviveImport() {
    gTest = "stored delays";
    hrtf::BuildOptions options;
    options.diffuseFieldEq = false;
    hrtf::RawHrirSet raw = rawImpulses(kRate, 256);
    raw.delays.resize(raw.directions.size() * 2);
    for (std::size_t d = 0; d < raw.directions.size(); ++d) {
        raw.delays[d * 2] = 3.0f;
        raw.delays[d * 2 + 1] = 13.0f;
    }
    std::unique_ptr<hrtf::HrtfSet> set = hrtf::buildSet(raw, options);
    CHECK(set != nullptr, "no set");
    if (!set) return;
    std::printf("  stored 10 samples between the ears, built %.2f\n", static_cast<double>(set->delays[1] - set->delays[0]));
    CHECK(std::fabs(set->delays[1] - set->delays[0] - 10.0f) < 0.01f, "%f", static_cast<double>(set->delays[1] - set->delays[0]));
    options.sampleRate = 96000.0f;
    set = hrtf::buildSet(raw, options);
    CHECK(set && std::fabs(set->delays[1] - set->delays[0] - 20.0f) < 0.02f, "not scaled with the rate");
}

// Finding of this round: the channel order came from ReceiverPosition alone, and CIPIC and FHK publish that with the opposite
// sign to their channels, so they imported with the ears swapped.
void testEarOrderComesFromTheData() {
    gTest = "ear order";
    constexpr std::uint32_t kTaps = 8;
    const float positions[4 * 3] = {0, 1, 0, 0, -1, 0, 1, 0.9f, 0, 1, -0.9f, 0};  // Left, right, and the same again further round.
    std::vector<float> leftFirst(4 * 2 * kTaps, 0.0f), rightFirst(4 * 2 * kTaps, 0.0f);
    for (int m = 0; m < 4; ++m) {
        const bool sourceLeft = positions[m * 3 + 1] > 0.0f;
        leftFirst[static_cast<std::size_t>(m) * 2 * kTaps] = sourceLeft ? 1.0f : 0.1f;
        leftFirst[(static_cast<std::size_t>(m) * 2 + 1) * kTaps] = sourceLeft ? 0.1f : 1.0f;
        rightFirst[static_cast<std::size_t>(m) * 2 * kTaps] = sourceLeft ? 0.1f : 1.0f;
        rightFirst[(static_cast<std::size_t>(m) * 2 + 1) * kTaps] = sourceLeft ? 1.0f : 0.1f;
    }
    for (const bool metadata : {false, true}) {
        CHECK(!hrtf::receiversSwapped(leftFirst.data(), positions, 4, kTaps, metadata), "the data say left first, metadata %d", metadata);
        CHECK(hrtf::receiversSwapped(rightFirst.data(), positions, 4, kTaps, metadata), "the data say right first, metadata %d", metadata);
    }
    const float frontal[3] = {1, 0, 0};
    CHECK(hrtf::receiversSwapped(leftFirst.data(), frontal, 1, kTaps, true), "no side sources: the metadata stands");
    CHECK(!hrtf::receiversSwapped(leftFirst.data(), frontal, 1, kTaps, false), "no side sources: the metadata stands");
}

// The same level fault as finding 13, on the import path: the build resamples before it does anything else.
void testImportRateKeepsLevel() {
    gTest = "import rate keeps level";
    hrtf::BuildOptions options;
    options.diffuseFieldEq = false;
    const auto dcOf = [](const hrtf::HrtfSet& set) {
        double sum = 0.0;
        for (std::uint32_t n = 0; n < set.taps; ++n) sum += set.filters[n];
        return sum;
    };
    const hrtf::RawHrirSet raw = rawImpulses(kRate, 256);
    const std::unique_ptr<hrtf::HrtfSet> native = hrtf::buildSet(raw, options);
    options.sampleRate = 96000.0f;
    const std::unique_ptr<hrtf::HrtfSet> doubled = hrtf::buildSet(raw, options);
    CHECK(native && doubled, "no set");
    if (!native || !doubled) return;
    const double difference = 20.0 * std::log10(std::fabs(dcOf(*doubled)) / std::fabs(dcOf(*native)));
    std::printf("  imported at 96 kHz against 48 kHz: %+.3f dB\n", difference);
    CHECK(std::fabs(difference) < 0.1, "%.2f dB", difference);
}

// Finding: resampling an impulse response as a waveform moved its level by outRate / inRate.
void testHrtfRateKeepsLevel() {
    gTest = "HRTF rate keeps level";
    const auto level = [](float rate, std::uint32_t voices) {
        ms_engine_desc d;
        ms_engine_desc_default(&d);
        d.sample_rate = rate;
        d.block_frames = rate > 50000.0f ? 256 : 128;
        d.max_hrtf_voices = voices;
        ms_engine* e = nullptr;
        OK(ms_engine_create(&d, &e));
        ms_hrtf* set = nullptr;
        OK(ms_hrtf_load(dataPath("mit_kemar.mshrtf").c_str(), rate, &set));
        OK(ms_engine_set_hrtf(e, set));
        const std::vector<float> tone = sineOf(static_cast<std::size_t>(rate), 1000.0, rate);
        ms_buffer b;
        ms_source s;
        OK(ms_buffer_create(e, tone.data(), tone.size(), 1, rate, &b));
        OK(ms_source_create(e, &s));
        const ms_source_params p = plainParams(horizontal(45.0, 1.2));
        OK(ms_source_set_params(e, s, &p));
        OK(ms_source_set_buffer(e, s, b, 0, 1));
        OK(ms_source_play(e, s, 0));
        OK(ms_engine_commit(e));
        const auto frames = static_cast<std::uint32_t>(rate / 2.0f);
        const Render r = renderFrames(e, frames);
        ms_hrtf_release(set);
        ms_engine_destroy(e);
        return rmsDb(r.l, r.r, frames / 2);
    };
    for (const std::uint32_t voices : {8u, 0u}) {
        const double reference = level(48000.0f, voices);
        for (const float rate : {44100.0f, 96000.0f}) {
            const double difference = level(rate, voices) - reference;
            std::printf("  1 kHz at %.0f Hz with %s against 48 kHz: %+.3f dB\n", static_cast<double>(rate), voices ? "direct voices" : "the ambisonic mix", difference);
            CHECK(std::fabs(difference) < 0.25, "%.2f dB at %f Hz", difference, static_cast<double>(rate));
        }
    }
}

void testDiffuseField() {
    gTest = "diffuse field";
    for (const auto& [name, set] : gMeasured) {
        const auto* s = reinterpret_cast<const hrtf::HrtfSet*>(set);
        std::vector<double> power(257, 0.0);
        core::RealFft fft(512);
        std::vector<float> time(512), re(256), im(256);
        for (std::uint32_t d = 0; d < s->directions; ++d)
            for (std::uint32_t ear = 0; ear < 2; ++ear) {
                std::fill(time.begin(), time.end(), 0.0f);
                const float* h = s->filters.data() + (d * 2ull + ear) * s->taps;
                std::copy(h, h + s->taps, time.begin());
                fft.forward(time.data(), re.data(), im.data());
                for (std::uint32_t k = 1; k < 256; ++k) power[k] += 0.5 * s->quadrature[d] * 0.25 * (static_cast<double>(re[k]) * re[k] + static_cast<double>(im[k]) * im[k]);
            }
        double low = 1e9, high = -1e9;
        for (std::uint32_t k = 3; k < 192; ++k) {  // 280 Hz to 18 kHz
            const double db = 10.0 * std::log10(power[k]);
            low = std::min(low, db);
            high = std::max(high, db);
        }
        std::printf("  %s: diffuse-field average spans %.2f to %.2f dB, 280 Hz to 18 kHz\n", name.c_str(), low, high);
        CHECK(high - low < 6.0, "%s spread %f dB", name.c_str(), high - low);
    }
}

void testConcurrentCommits() {
    gTest = "lock-free boundary";
    // Snapshots are whole: a writer publishes pairs that must always match, and the reader takes them while the writer runs.
    core::TripleIndex index;
    static std::uint64_t slots[3][2];
    std::atomic<bool> stop{false};
    std::thread writer([&] {
        std::uint64_t value = 1;
        while (!stop.load(std::memory_order_relaxed)) {
            const std::uint32_t w = index.writeIndex();
            slots[w][0] = value;
            std::atomic_thread_fence(std::memory_order_relaxed);
            slots[w][1] = value;
            index.publish();
            ++value;
        }
    });
    std::uint64_t torn = 0, last = 0, regressions = 0, fresh = 0;
    for (int i = 0; i < 20000000 && fresh < 200000; ++i) {
        if (!index.acquire()) continue;
        ++fresh;
        const std::uint32_t r = index.readIndex();
        const std::uint64_t a = slots[r][0], b = slots[r][1];
        torn += a != b;
        regressions += a < last;
        last = a;
    }
    stop = true;
    writer.join();
    std::printf("  %llu snapshots read while written: %llu torn, %llu out of order\n", static_cast<unsigned long long>(fresh), static_cast<unsigned long long>(torn), static_cast<unsigned long long>(regressions));
    CHECK(torn == 0 && regressions == 0, "torn %llu regress %llu", static_cast<unsigned long long>(torn), static_cast<unsigned long long>(regressions));
    CHECK(fresh > 1000, "only %llu snapshots", static_cast<unsigned long long>(fresh));

    core::SpscRing<std::uint64_t> ring;
    ring.init(64);
    std::atomic<bool> done{false};
    std::thread producer([&] {
        for (std::uint64_t i = 1; i <= 2000000;) {
            if (ring.push(i)) ++i;
        }
        done = true;
    });
    std::uint64_t expected = 1, misordered = 0, value;
    while (expected <= 2000000) {
        if (ring.pop(value)) {
            misordered += value != expected;
            ++expected;
        }
    }
    producer.join();
    CHECK(misordered == 0, "%llu misordered", static_cast<unsigned long long>(misordered));
}

// Flops of one real FFT of n points: Accelerate's split-radix class, 2.5 n log2 n is a standard upper count.
double fftFlops(std::uint32_t n) { return 2.5 * n * std::log2(static_cast<double>(n)); }

void testBudget() {
    gTest = "budget";
    // M3 performance core: 4.05 GHz, four 128-bit floating point pipes, so 64.8 GFLOP/s single precision at the absolute peak.
    // Real sustained DSP reaches a fraction of that; 15% is used here. Budget: half of one performance core.
    constexpr double kAssumedFlopsPerSecond = 64.8e9 * 0.15, kBudgetShare = 0.5;
    struct Case {
        const char* label;
        std::uint32_t direct, shared, order, reflections;
        bool moving;
    };
    const Case cases[] = {
        {"32 direct, moving, order 3 mix", 32, 0, 3, 0, true},
        {"64 direct, moving, order 3 mix", 64, 0, 3, 0, true},
        {"64 direct + 192 shared, moving", 64, 192, 3, 0, true},
        {"32 direct with 12 reflections each", 32, 0, 3, 12, true},
        {"64 direct, still, order 3 mix", 64, 0, 3, 0, false},
    };
    std::vector<float> noise(48000);
    std::mt19937 rng(13);
    for (float& v : noise) v = std::uniform_real_distribution<float>(-0.3f, 0.3f)(rng);
    std::printf("  counted operations per block of 128 at 48 kHz (block budget %.2f ms, half core %.0f kflop assumed)\n", 128.0 / 48.0, kBudgetShare * kAssumedFlopsPerSecond * 128.0 / 48000.0 / 1e3);
    std::printf("    %-38s %7s %7s %9s %9s %8s %8s\n", "case", "fft", "ifft", "cmac", "mac", "kflop", "M3 %%core");
    for (const Case& c : cases) {
        const std::uint32_t total = c.direct + c.shared;
        ms_engine_desc d;
        ms_engine_desc_default(&d);
        d.max_sources = std::max(total, 1u);
        d.max_hrtf_voices = c.direct;
        d.ambisonic_order = c.order;
        ms_engine* e = nullptr;
        OK(ms_engine_create(&d, &e));
        OK(ms_engine_set_hrtf(e, gSphere));
        ms_buffer b;
        OK(ms_buffer_create(e, noise.data(), noise.size(), 1, kRate, &b));
        std::vector<ms_source> ids(total);
        for (std::uint32_t i = 0; i < total; ++i) {
            OK(ms_source_create(e, &ids[i]));
            ms_source_params p;
            ms_source_params_default(&p);
            p.position = horizontal(i * 37.0, 2.0f + static_cast<float>(i % 7));
            p.velocity = {3.0f, 0.0f, 1.0f};
            p.lowpass_hz = 8000.0f;
            OK(ms_source_set_params(e, ids[i], &p));
            OK(ms_source_set_buffer(e, ids[i], b, 0, 1));
            OK(ms_source_play(e, ids[i], 0));
            if (c.reflections) {
                ms_reflection_tap taps[MS_MAX_REFLECTION_TAPS];
                for (std::uint32_t t = 0; t < c.reflections; ++t) taps[t] = {0.005f + 0.004f * static_cast<float>(t), horizontal(t * 53.0, 1.0f), {0.4f, 0.3f, 0.2f}};
                OK(ms_source_set_reflections(e, ids[i], taps, c.reflections));
            }
        }
        OK(ms_engine_commit(e));
        // Settle, then measure across a span with a new pose every 60 Hz frame, the worst case for filter resynthesis.
        (void)renderFrames(e, 2560);
        ms_op_counts before, after;
        ms_engine_get_op_counts(e, &before);
        double wall = 0.0;
        std::vector<float> l(128), r(128);
        constexpr std::uint32_t kBlocks = 2000;
        for (std::uint32_t block = 0; block < kBlocks; ++block) {
            if (c.moving && block % 6 == 0) {
                for (std::uint32_t i = 0; i < total; ++i)
                    OK(ms_source_set_pose(e, ids[i], horizontal(i * 37.0 + block * 0.5, 2.0f + static_cast<float>(i % 7)), {3.0f, 0.0f, 1.0f}));
                OK(ms_engine_commit(e));
            }
            const auto start = std::chrono::steady_clock::now();
            ms_engine_render(e, l.data(), r.data(), 128);
            wall += std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
        }
        ms_engine_get_op_counts(e, &after);
        const double fft = static_cast<double>(after.fft_forward - before.fft_forward) / kBlocks;
        const double ifft = static_cast<double>(after.fft_inverse - before.fft_inverse) / kBlocks;
        const double cmac = static_cast<double>(after.spectrum_macs - before.spectrum_macs) / kBlocks;
        const double mac = static_cast<double>(after.real_macs - before.real_macs) / kBlocks;
        const double flops = (fft + ifft) * fftFlops(256) + cmac * 8.0 + mac * 2.0;
        const double share = flops / (kAssumedFlopsPerSecond * 128.0 / 48000.0);
        std::printf("    %-38s %7.0f %7.0f %9.0f %9.0f %8.0f %8.1f   (this host, wall: %.1f%% of realtime)\n", c.label, fft, ifft, cmac, mac, flops / 1e3, share * 100.0,
                    wall / kBlocks / (128.0 / 48000.0) * 100.0);
        CHECK(share < kBudgetShare, "%s over budget: %f of a core", c.label, share);
        ms_engine_destroy(e);
    }
}

// Energy per channel of a speaker render of one source or bed.
std::vector<double> speakerEnergy(const ms_speaker_layout& layout, ms_vec3 position, bool bed, ms_quat head, std::uint32_t frames,
                                  std::vector<std::vector<float>>* keep = nullptr, std::uint32_t chunk = 0) {
    ms_engine_desc d;
    ms_engine_desc_default(&d);
    d.smoothing_seconds = 0.0f;
    d.speakers = layout;
    ms_engine* e = nullptr;
    if (ms_engine_create(&d, &e) != MS_OK) std::abort();
    std::vector<float> noise(frames * 4);
    std::mt19937 rng(11);
    std::normal_distribution<float> normal(0.0f, 0.1f);
    ms_buffer b;
    if (bed) {
        // A first-order plane wave from position, as a direction in listener axes: AmbiX X front, Y left, Z up.
        const float n = std::sqrt(position.x * position.x + position.y * position.y + position.z * position.z);
        const float x = -position.z / n, y = -position.x / n, z = position.y / n;
        for (std::uint32_t t = 0; t < frames; ++t) {
            const float s = normal(rng);
            noise[4 * t] = s;
            noise[4 * t + 1] = s * y;
            noise[4 * t + 2] = s * z;
            noise[4 * t + 3] = s * x;
        }
        OK(ms_buffer_create(e, noise.data(), frames, 4, kRate, &b));
        ms_bed bd;
        OK(ms_bed_create(e, &bd));
        OK(ms_bed_set_buffer(e, bd, b, 1));
        OK(ms_bed_play(e, bd, 0));
    } else {
        for (std::uint32_t t = 0; t < frames; ++t) noise[t] = normal(rng);
        OK(ms_buffer_create(e, noise.data(), frames, 1, kRate, &b));
        ms_source s;
        OK(ms_source_create(e, &s));
        const ms_source_params p = plainParams(position);
        OK(ms_source_set_params(e, s, &p));
        OK(ms_source_set_buffer(e, s, b, 0, 1));
        OK(ms_source_play(e, s, 0));
    }
    ms_head_set_orientation(e, head);
    OK(ms_engine_commit(e));
    const std::uint32_t n = ms_engine_output_channels(e);
    std::vector<std::vector<float>> out(n, std::vector<float>(frames));
    std::vector<float*> pointers(n);
    for (std::uint32_t done = 0; done < frames;) {
        const std::uint32_t take = chunk ? std::min(chunk, frames - done) : frames;
        for (std::uint32_t c = 0; c < n; ++c) pointers[c] = out[c].data() + done;
        ms_engine_render_speakers(e, pointers.data(), take);
        done += take;
    }
    ms_engine_destroy(e);
    std::vector<double> energy(n, 0.0);
    for (std::uint32_t c = 0; c < n; ++c)
        for (std::uint32_t t = 1024; t < frames; ++t) energy[c] += static_cast<double>(out[c][t]) * out[c][t];
    if (keep) *keep = std::move(out);
    return energy;
}

void testSpeakers() {
    gTest = "speakers";
    const std::pair<ms_speaker_preset, std::uint32_t> presets[] = {{MS_SPEAKERS_MONO, 1}, {MS_SPEAKERS_STEREO, 2}, {MS_SPEAKERS_QUAD, 4},
                                                                   {MS_SPEAKERS_5_1, 6},  {MS_SPEAKERS_7_1, 8},     {MS_SPEAKERS_5_1_4, 10},
                                                                   {MS_SPEAKERS_7_1_4, 12}};
    std::mt19937 rng(5);
    std::normal_distribution<float> normal;
    for (const auto& [preset, count] : presets) {
        ms_speaker_layout layout;
        OK(ms_speaker_layout_preset(preset, &layout));
        CHECK(layout.count == count, "preset %d has %u speakers", preset, layout.count);
        core::SpeakerLayout built;
        CHECK(built.build(layout, 3), "preset %d builds", preset);
        float gains[core::kMaxSpeakers];
        // A source at a speaker plays from that speaker alone.
        double worstAt = 0.0;
        for (std::uint32_t s = 0; s < layout.count; ++s) {
            if (layout.speakers[s].lfe) continue;
            built.pan(core::speakerDirection(layout.speakers[s].azimuth_degrees, layout.speakers[s].elevation_degrees), gains);
            for (std::uint32_t k = 0; k < layout.count; ++k) worstAt = std::max(worstAt, std::fabs(static_cast<double>(gains[k]) - (k == s ? 1.0 : 0.0)));
        }
        CHECK(worstAt < 1e-4, "preset %d: a source at a speaker spills %.2e", preset, worstAt);
        // Energy is one from everywhere, and the LFE takes none.
        double worstEnergy = 0.0, worstAngle = 0.0, lfe = 0.0;
        for (int trial = 0; trial < 2000; ++trial) {
            const core::Vec3 d = core::normalized(core::Vec3{normal(rng), normal(rng), normal(rng)});
            built.pan(d, gains);
            double energy = 0.0;
            for (std::uint32_t k = 0; k < layout.count; ++k) {
                energy += gains[k] * gains[k];
                if (layout.speakers[k].lfe) lfe += std::fabs(gains[k]);
            }
            worstEnergy = std::max(worstEnergy, std::fabs(energy - 1.0));
        }
        // On a ring, the panned image (the gain-weighted sum of speaker directions) points where the source is, wherever two neighbours
        // hold it between them: everywhere on the horizon but behind a stereo pair.
        if (preset != MS_SPEAKERS_MONO && preset != MS_SPEAKERS_5_1_4 && preset != MS_SPEAKERS_7_1_4) {
            for (int a = -180; a < 180; ++a) {
                if (preset == MS_SPEAKERS_STEREO && std::abs(a) > 30) continue;
                const core::Vec3 d = core::speakerDirection(static_cast<float>(a), 0.0f);
                built.pan(d, gains);
                core::Vec3 image{};
                for (std::uint32_t k = 0; k < layout.count; ++k)
                    if (!layout.speakers[k].lfe)
                        image = image + core::speakerDirection(layout.speakers[k].azimuth_degrees, layout.speakers[k].elevation_degrees) * gains[k];
                worstAngle = std::max(worstAngle, std::acos(std::clamp(static_cast<double>(core::dot(core::normalized(image), d)), -1.0, 1.0)) * 180.0 / M_PI);
            }
        }
        std::printf("  %2u speakers: energy within %.1e of one, image within %.2f degrees on the horizon\n", layout.count, worstEnergy, worstAngle);
        CHECK(worstEnergy < 1e-4 && lfe == 0.0 && worstAngle < 0.5, "preset %d pans", preset);

        // AllRAD at order 3: a plane wave from anywhere carries about the same energy, and the image (energy vector) points the right way.
        float encode[core::kMaxChannels];
        // 5.1 leaves 140 degrees open behind; its image is measured over the front half, where its speakers are.
        const bool irregular = preset == MS_SPEAKERS_5_1 || preset == MS_SPEAKERS_5_1_4;
        double lowest = INFINITY, highest = 0.0, worstImage = 0.0;
        for (int trial = 0; trial < 400; ++trial) {
            const float az = static_cast<float>(trial % 40) * 9.0f, el = preset == MS_SPEAKERS_5_1_4 || preset == MS_SPEAKERS_7_1_4 ? static_cast<float>(trial / 40) * 5.0f : 0.0f;
            const core::Vec3 d = core::speakerDirection(az, el);
            core::encodeDirection(core::headToAmbisonic(d), 3, encode);
            double energy = 0.0;
            core::Vec3 image{};
            for (std::uint32_t k = 0; k < layout.count; ++k) {
                const float* row = built.decodeRow(k);
                double feed = 0.0;
                for (std::uint32_t c = 0; c < built.decodeChannels(); ++c) feed += row[c] * encode[c];
                energy += feed * feed;
                if (!layout.speakers[k].lfe)
                    image = image + core::speakerDirection(layout.speakers[k].azimuth_degrees, layout.speakers[k].elevation_degrees) * static_cast<float>(feed * feed);
            }
            lowest = std::min(lowest, energy);
            highest = std::max(highest, energy);
            if (preset != MS_SPEAKERS_MONO && preset != MS_SPEAKERS_STEREO && (!irregular || az <= 90.0f || az >= 270.0f))
                worstImage = std::max(worstImage, std::acos(std::clamp(static_cast<double>(core::dot(core::normalized(image), d)), -1.0, 1.0)) * 180.0 / M_PI);
        }
        const double spread = 10.0 * std::log10(highest / lowest);
        std::printf("     decoder: plane-wave energy spread %.2f dB, energy image within %.1f degrees%s\n", spread, worstImage, irregular ? " in front" : "");
        if (preset != MS_SPEAKERS_MONO && preset != MS_SPEAKERS_STEREO) CHECK(spread < 4.0 && worstImage < 15.0, "preset %d decodes", preset);
    }
    ms_speaker_layout bad{};
    bad.count = 2;
    bad.speakers[0] = {0, 0, 1};
    bad.speakers[1] = {0, 0, 1};
    float gains[core::kMaxSpeakers];
    CHECK(ms_speaker_layout_pan(&bad, {0, 0, -1}, gains) == MS_ERR_INVALID_ARGUMENT, "a layout of LFE alone is refused");

    // The engine on 5.1: a source in front plays from the centre, one at 110 degrees left from the left surround, and the LFE stays silent.
    ms_speaker_layout surround;
    OK(ms_speaker_layout_preset(MS_SPEAKERS_5_1, &surround));
    const auto loudest = [](const std::vector<double>& energy) { return static_cast<std::uint32_t>(std::max_element(energy.begin(), energy.end()) - energy.begin()); };
    const std::vector<double> front = speakerEnergy(surround, {0, 0, -2}, false, {0, 0, 0, 1}, 9600);
    CHECK(loudest(front) == 2 && front[3] == 0.0, "front source: loudest %u, LFE %g", loudest(front), front[3]);
    double others = 0.0;
    for (std::uint32_t k = 0; k < front.size(); ++k) others += k == 2 ? 0.0 : front[k];
    CHECK(others < 1e-6 * front[2], "a source on the centre speaker leaks %.2e into the others", others / front[2]);
    const core::Vec3 left = core::speakerDirection(110.0f, 0.0f) * 2.0f;
    const std::vector<double> surroundLeft = speakerEnergy(surround, {left.x, left.y, left.z}, false, {0, 0, 0, 1}, 9600);
    CHECK(loudest(surroundLeft) == 4, "a source at 110 degrees left plays from the left surround, not %u", loudest(surroundLeft));
    // Turning the head changes nothing: the speakers stand in the room.
    const float h = std::sqrt(0.5f);
    const std::vector<double> turned = speakerEnergy(surround, {0, 0, -2}, false, {0, h, 0, h}, 9600);
    double change = 0.0;
    for (std::uint32_t k = 0; k < front.size(); ++k) change = std::max(change, std::fabs(turned[k] - front[k]));
    CHECK(change == 0.0, "head orientation moved the speaker image by %g", change);
    // A first-order bed from the left decodes onto the left speakers.
    const std::vector<double> bedLeft = speakerEnergy(surround, {-1, 0, 0}, true, {0, 0, 0, 1}, 9600);
    CHECK(bedLeft[0] > bedLeft[1] && bedLeft[4] > bedLeft[5] && bedLeft[3] == 0.0, "bed from the left: L %g R %g Ls %g Rs %g", bedLeft[0], bedLeft[1], bedLeft[4], bedLeft[5]);
    // Chunking does not change a sample.
    std::vector<std::vector<float>> whole, pieces;
    (void)speakerEnergy(surround, {1, 0.5f, -2}, false, {0, 0, 0, 1}, 9600, &whole);
    (void)speakerEnergy(surround, {1, 0.5f, -2}, false, {0, 0, 0, 1}, 9600, &pieces, 77);
    CHECK(whole == pieces, "speaker output depends on how frames are chunked");
    // The stereo entry point on a speaker engine folds the same render down.
    ms_engine_desc d;
    ms_engine_desc_default(&d);
    d.speakers = surround;
    ms_engine* e = nullptr;
    OK(ms_engine_create(&d, &e));
    CHECK(ms_engine_output_channels(e) == 6, "5.1 renders six channels");
    std::vector<float> noise(4800);
    for (float& v : noise) v = 0.1f * normal(rng);
    ms_buffer b;
    ms_source s;
    OK(ms_buffer_create(e, noise.data(), noise.size(), 1, kRate, &b));
    OK(ms_source_create(e, &s));
    const ms_source_params p = plainParams({-2, 0, 0});
    OK(ms_source_set_params(e, s, &p));
    OK(ms_source_set_buffer(e, s, b, 0, 1));
    OK(ms_source_play(e, s, 0));
    OK(ms_engine_commit(e));
    Render folded = renderFrames(e, 4800);
    double fl = 0.0, fr = 0.0;
    for (std::size_t t = 1024; t < folded.l.size(); ++t) {
        fl += folded.l[t] * folded.l[t];
        fr += folded.r[t] * folded.r[t];
    }
    CHECK(fl > 100.0 * fr && fl > 0.0, "a source on the left folds down to the left: %g against %g", fl, fr);
    ms_engine_destroy(e);
    ms_engine_desc_default(&d);
    d.speakers.count = 2;
    d.speakers.speakers[0] = {0, 0, 1};
    d.speakers.speakers[1] = {10, 0, 1};
    CHECK(ms_engine_create(&d, &e) == MS_ERR_INVALID_ARGUMENT, "a layout of LFE channels alone is refused");
}

// A channel layout with descriptions, as a device would report one.
std::vector<unsigned char> describedLayout(const std::vector<AudioChannelDescription>& descriptions) {
    std::vector<unsigned char> bytes(offsetof(AudioChannelLayout, mChannelDescriptions) + descriptions.size() * sizeof(AudioChannelDescription));
    auto* layout = reinterpret_cast<AudioChannelLayout*>(bytes.data());
    layout->mChannelLayoutTag = kAudioChannelLayoutTag_UseChannelDescriptions;
    layout->mNumberChannelDescriptions = static_cast<UInt32>(descriptions.size());
    std::memcpy(layout->mChannelDescriptions, descriptions.data(), descriptions.size() * sizeof(AudioChannelDescription));
    return bytes;
}

void testDeviceChannels() {
    gTest = "device channels";
    const auto near = [](float a, float b) { return std::fabs(a - b) < 1e-3f; };
    // A layout named by its tag: 5.1 in the MPEG order.
    AudioChannelLayout tagged{};
    tagged.mChannelLayoutTag = kAudioChannelLayoutTag_MPEG_5_1_A;
    ms_speaker_layout layout;
    CHECK(core::speakerLayoutFromChannelLayout(tagged, 6, layout) && layout.count == 6, "5.1 by tag");
    const float azimuth51[] = {30, -30, 0, 0, 110, -110};
    for (std::uint32_t c = 0; c < 6; ++c)
        CHECK(c == 3 ? layout.speakers[c].lfe != 0 : (!layout.speakers[c].lfe && near(layout.speakers[c].azimuth_degrees, azimuth51[c])), "5.1 channel %u at %g", c,
              layout.speakers[c].azimuth_degrees);
    core::SpeakerLayout built;
    CHECK(built.build(layout, 3), "the device's 5.1 builds");
    // Dolby Atmos 7.1.4 by tag: twelve channels, four of them up.
    tagged.mChannelLayoutTag = kAudioChannelLayoutTag_Atmos_7_1_4;
    CHECK(core::speakerLayoutFromChannelLayout(tagged, 12, layout) && layout.count == 12, "7.1.4 by tag");
    std::uint32_t up = 0, placed = 0;
    for (std::uint32_t c = 0; c < 12; ++c) {
        placed += !layout.speakers[c].lfe;
        up += !layout.speakers[c].lfe && layout.speakers[c].elevation_degrees > 0.0f;
    }
    CHECK(placed == 11 && up == 4 && built.build(layout, 3), "7.1.4: %u placed, %u up", placed, up);
    // Descriptions: coordinates (Core Audio's azimuth is positive to the right), a repeated position, a discrete channel, and a device with
    // more channels than the layout describes.
    AudioChannelDescription right45{kAudioChannelLabel_UseCoordinates, kAudioChannelFlags_SphericalCoordinates, {45, 20, 1}};
    AudioChannelDescription left{kAudioChannelLabel_Left, 0, {}}, again{kAudioChannelLabel_LeftTotal, 0, {}}, discrete{kAudioChannelLabel_Discrete_0, 0, {}};
    AudioChannelDescription box{kAudioChannelLabel_UseCoordinates, kAudioChannelFlags_RectangularCoordinates, {-1, 1, 0}};
    const auto bytes = describedLayout({right45, left, again, discrete, box});
    CHECK(core::speakerLayoutFromChannelLayout(*reinterpret_cast<const AudioChannelLayout*>(bytes.data()), 7, layout) && layout.count == 7, "described layout");
    CHECK(near(layout.speakers[0].azimuth_degrees, -45) && near(layout.speakers[0].elevation_degrees, 20) && !layout.speakers[0].lfe, "spherical coordinates");
    CHECK(near(layout.speakers[1].azimuth_degrees, 30) && layout.speakers[2].lfe && layout.speakers[3].lfe, "a repeated place and a discrete channel stay silent");
    CHECK(near(layout.speakers[4].azimuth_degrees, 45) && near(layout.speakers[4].elevation_degrees, 0), "rectangular coordinates: front left");
    CHECK(layout.speakers[5].lfe && layout.speakers[6].lfe, "channels past the layout stay silent");
    const auto nothing = describedLayout({discrete});
    CHECK(!core::speakerLayoutFromChannelLayout(*reinterpret_cast<const AudioChannelLayout*>(nothing.data()), 1, layout), "a layout that places nothing");

    // Writing: one interleaved buffer of eight channels with six to fill, then one buffer per channel.
    constexpr std::uint32_t frames = 16;
    std::vector<std::vector<float>> source(6, std::vector<float>(frames));
    std::vector<const float*> pointers;
    for (std::uint32_t c = 0; c < 6; ++c) {
        for (std::uint32_t t = 0; t < frames; ++t) source[c][t] = static_cast<float>(c * 100 + t);
        pointers.push_back(source[c].data());
    }
    std::vector<float> interleaved(8 * frames, -1.0f);
    AudioBufferList one{1, {{8, static_cast<UInt32>(interleaved.size() * sizeof(float)), interleaved.data()}}};
    core::writeDeviceChannels(pointers.data(), 6, frames, one);
    bool right = true;
    for (std::uint32_t t = 0; t < frames; ++t)
        for (std::uint32_t k = 0; k < 8; ++k) right = right && interleaved[t * 8 + k] == (k < 6 ? static_cast<float>(k * 100 + t) : 0.0f);
    CHECK(right, "interleaved channels land in place, the rest silent");
    std::vector<unsigned char> listBytes(offsetof(AudioBufferList, mBuffers) + 7 * sizeof(AudioBuffer));
    auto* split = reinterpret_cast<AudioBufferList*>(listBytes.data());
    split->mNumberBuffers = 7;
    std::vector<std::vector<float>> planes(7, std::vector<float>(frames, -1.0f));
    for (std::uint32_t b = 0; b < 7; ++b) split->mBuffers[b] = {1, frames * static_cast<UInt32>(sizeof(float)), planes[b].data()};
    core::writeDeviceChannels(pointers.data(), 6, frames, *split);
    right = true;
    for (std::uint32_t b = 0; b < 7; ++b)
        for (std::uint32_t t = 0; t < frames; ++t) right = right && planes[b][t] == (b < 6 ? static_cast<float>(b * 100 + t) : 0.0f);
    CHECK(right, "a buffer per channel: each its own, the extra silent");
}

std::complex<double> tapResponse(const std::vector<float>& taps, double hz) {
    std::complex<double> sum = 0.0;
    for (std::size_t n = 0; n < taps.size(); ++n) sum += static_cast<double>(taps[n]) * std::polar(1.0, -2.0 * M_PI * hz * static_cast<double>(n) / kRate);
    return sum;
}

// Renders a noise source at a position through a fresh engine, with the canceller on from the start when params is given, and switched
// off after offFrames when that is nonzero.
Render crosstalkRender(ms_vec3 position, const ms_crosstalk_params* params, std::uint32_t frames, std::uint32_t offFrames = 0) {
    ms_engine* e = makeEngine(gSphere);
    std::vector<float> noise(frames);
    std::mt19937 rng(17);
    std::normal_distribution<float> normal(0.0f, 0.1f);
    for (float& v : noise) v = normal(rng);
    ms_buffer b;
    ms_source s;
    OK(ms_buffer_create(e, noise.data(), noise.size(), 1, kRate, &b));
    OK(ms_source_create(e, &s));
    const ms_source_params p = plainParams(position);
    OK(ms_source_set_params(e, s, &p));
    OK(ms_source_set_buffer(e, s, b, 0, 0));
    OK(ms_source_play(e, s, 0));
    if (params) OK(ms_engine_set_crosstalk(e, params));
    OK(ms_engine_commit(e));
    Render out;
    out.l.resize(frames);
    out.r.resize(frames);
    if (offFrames) {
        ms_engine_render(e, out.l.data(), out.r.data(), offFrames);
        OK(ms_engine_set_crosstalk(e, nullptr));
        CHECK(ms_engine_crosstalk_latency_frames(e) == 0, "no latency once off");
        OK(ms_engine_commit(e));
        ms_engine_render(e, out.l.data() + offFrames, out.r.data() + offFrames, frames - offFrames);
    } else {
        if (params) CHECK(ms_engine_crosstalk_latency_frames(e) == 256, "latency %u", ms_engine_crosstalk_latency_frames(e));
        ms_engine_render(e, out.l.data(), out.r.data(), frames);
    }
    ms_engine_destroy(e);
    return out;
}

void testCrosstalk() {
    gTest = "crosstalk";
    ms_crosstalk_params params;
    ms_crosstalk_params_default(&params);
    for (float span : {10.0f, 15.0f, 30.0f}) {
        params.speaker_span_degrees = span;
        const auto filter = core::designCrosstalk(params, kRate, 128, kSpeed);
        CHECK(filter != nullptr, "span %g designs", span);
        if (!filter) continue;
        // In the band: the far ear's share of one channel, with and without the canceller, through the same spherical-head plant.
        double rejection = 0.0, plain = 0.0, worstBoost = -INFINITY;
        int count = 0;
        for (double hz = 300.0; hz <= 5000.0; hz *= std::pow(2.0, 1.0 / 6.0), ++count) {
            const std::complex<double> r = core::crosstalkPlant(params, hz, kSpeed);
            const std::complex<double> cs = tapResponse(filter->sum, hz), cd = tapResponse(filter->difference, hz);
            const std::complex<double> near = ((1.0 + r) * cs + (1.0 - r) * cd) / 2.0, far = ((1.0 + r) * cs - (1.0 - r) * cd) / 2.0;
            rejection += 20.0 * std::log10(std::abs(near) / std::abs(far));
            plain += 20.0 * std::log10(1.0 / std::abs(r));
        }
        rejection /= count;
        plain /= count;
        for (double hz = 20.0; hz <= 20000.0; hz *= std::pow(2.0, 1.0 / 12.0)) {
            const double cs = std::abs(tapResponse(filter->sum, hz)), cd = std::abs(tapResponse(filter->difference, hz));
            worstBoost = std::max(worstBoost, 20.0 * std::log10(std::max(cs, cd)));
        }
        double outside = 0.0;
        for (double hz : {40.0, 60.0, 80.0, 15000.0, 18000.0})
            for (const auto* taps : {&filter->sum, &filter->difference}) {
                const double db = 20.0 * std::log10(std::abs(tapResponse(*taps, hz)));
                outside = std::max(outside, std::fabs(db));
            }
        std::printf("  span %4.1f: far ear %5.1f dB down with the canceller, %4.1f dB without; largest boost %4.1f dB; out of band within %.2f dB\n", span,
                    rejection, plain, worstBoost, outside);
        CHECK(rejection > 20.0 && rejection - plain > 12.0, "span %g: %g dB against %g", span, rejection, plain);
        CHECK(worstBoost < params.max_boost_db + 0.5, "span %g boosts %g dB", span, worstBoost);
        CHECK(outside < 1.0, "span %g: out of band %g dB from plain stereo", span, outside);
    }
    // Off the design: a canceller made for 15 degrees and an average head, heard with the speakers at 12 or 18 degrees (the player leaning
    // in or back) or by a smaller head. It should still beat plain stereo, though by less.
    {
        ms_crosstalk_params designed;
        ms_crosstalk_params_default(&designed);
        const auto filter = core::designCrosstalk(designed, kRate, 128, kSpeed);
        for (const auto& [span, radius] : {std::pair{12.0f, 0.0875f}, std::pair{18.0f, 0.0875f}, std::pair{15.0f, 0.075f}}) {
            ms_crosstalk_params actual = designed;
            actual.speaker_span_degrees = span;
            actual.head_radius_m = radius;
            double rejection = 0.0, plain = 0.0;
            int count = 0;
            for (double hz = 300.0; hz <= 5000.0; hz *= std::pow(2.0, 1.0 / 6.0), ++count) {
                const std::complex<double> r = core::crosstalkPlant(actual, hz, kSpeed);
                const std::complex<double> cs = tapResponse(filter->sum, hz), cd = tapResponse(filter->difference, hz);
                rejection += 20.0 * std::log10(std::abs((1.0 + r) * cs + (1.0 - r) * cd) / std::abs((1.0 + r) * cs - (1.0 - r) * cd));
                plain += 20.0 * std::log10(1.0 / std::abs(r));
            }
            std::printf("  designed for 15 degrees, heard at %4.1f with a %.1f cm head: far ear %5.1f dB down, %4.1f dB plain\n", span, radius * 100.0,
                        rejection / count, plain / count);
            CHECK(rejection / count - plain / count > 6.0, "off the design at %g degrees: %g against %g", span, rejection / count, plain / count);
        }
    }
    params.speaker_span_degrees = 70.0f;
    CHECK(core::designCrosstalk(params, kRate, 128, kSpeed) == nullptr, "a span past 60 degrees is refused");
    ms_engine* bad = makeEngine(gSphere);
    CHECK(ms_engine_set_crosstalk(bad, &params) == MS_ERR_INVALID_ARGUMENT, "and so through the API");
    ms_engine_destroy(bad);

    // The render path: the engine's output equals the plain render shuffled, filtered and recombined offline.
    ms_crosstalk_params_default(&params);
    const std::uint32_t frames = 9600;
    const Render dry = crosstalkRender({-2, 0, -1}, nullptr, frames), wet = crosstalkRender({-2, 0, -1}, &params, frames);
    const auto filter = core::designCrosstalk(params, kRate, 128, kSpeed);
    double worst = 0.0, peak = 0.0;
    // The first block crosses over from plain stereo, as any change of filter does.
    for (std::uint32_t t = 128; t < frames; ++t) {
        double s = 0.0, d = 0.0;
        for (std::uint32_t k = 0; k < filter->taps && k <= t; ++k) {
            s += filter->sum[k] * (static_cast<double>(dry.l[t - k]) + dry.r[t - k]);
            d += filter->difference[k] * (static_cast<double>(dry.l[t - k]) - dry.r[t - k]);
        }
        worst = std::max({worst, std::fabs(wet.l[t] - 0.5 * (s + d)), std::fabs(wet.r[t] - 0.5 * (s - d))});
        peak = std::max(peak, static_cast<double>(std::fabs(wet.l[t])));
    }
    std::printf("  render path against offline filtering: worst error %.2e of a %.2e peak\n", worst, peak);
    CHECK(peak > 0.0 && worst < 1e-4 * peak + 1e-7, "render path differs by %g", worst);
    // Switched off, the master is exactly the plain render again from the next block on.
    const Render off = crosstalkRender({-2, 0, -1}, &params, frames, 4800);
    bool same = true;
    for (std::uint32_t t = 4992; t < frames; ++t) same = same && off.l[t] == dry.l[t] && off.r[t] == dry.r[t];
    CHECK(same, "switching the canceller off returns the plain render");
    // Speaker layouts have no use for it.
    ms_speaker_layout layout;
    OK(ms_speaker_layout_preset(MS_SPEAKERS_5_1, &layout));
    ms_engine_desc d;
    ms_engine_desc_default(&d);
    d.speakers = layout;
    ms_engine* e = nullptr;
    OK(ms_engine_create(&d, &e));
    ms_crosstalk_params_default(&params);
    CHECK(ms_engine_set_crosstalk(e, &params) == MS_ERR_UNSUPPORTED, "a speaker layout refuses the canceller");
    ms_engine_destroy(e);
}

}  // namespace

int main(int argumentCount, char** arguments) {
    std::setvbuf(stdout, nullptr, _IOLBF, 0);
    if (argumentCount > 1 && std::strcmp(arguments[1], "--sofa") == 0) {
        const char* sofa = std::getenv("MS_TEST_SOFA");
        if (!sofa) {
            std::printf("skipped: set MS_TEST_SOFA to a real file; none ships in the library\n");
            return 77;
        }
        if (ms_hrtf_create_sphere(kRadius, kRate, &gSphere) != MS_OK || !gSphere) return 2;
        testSofaFile(sofa);
        ms_hrtf_release(gSphere);
        std::printf("%d checks, %d failures\n", gChecks, gFailures);
        return gFailures ? 1 : 0;
    }
    const auto run = [](const char* name, const std::function<void()>& test) {
        const auto start = std::chrono::steady_clock::now();
        const int before = gFailures;
        std::printf("== %s\n", name);
        test();
        std::printf("   %s, %.2f s\n", gFailures == before ? "pass" : "FAIL", std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count());
    };
    testSphereModel();
    if (ms_hrtf_create_sphere(kRadius, kRate, &gSphere) != MS_OK || !gSphere) return 2;
    ms_hrtf_info info;
    ms_hrtf_get_info(gSphere, &info);
    std::printf("sphere set: %u directions, %u taps\n", info.directions, info.taps);
    for (const char* name : {"sadie2_d1_ku100.mshrtf", "sadie2_d2_kemar.mshrtf", "sadie2_h3.mshrtf", "mit_kemar.mshrtf"}) {
        ms_hrtf* set = nullptr;
        if (ms_hrtf_load(dataPath(name).c_str(), kRate, &set) == MS_OK) gMeasured.emplace_back(name, set);
        else std::printf("not loaded: %s\n", name);
    }
    CHECK(gMeasured.size() >= 2, "vendored sets missing: %zu", gMeasured.size());
    run("levels", testFilterIdentityAndLevels);
    run("head model, analytic sphere", [] { testHeadModelCues("sphere", gSphere, true); });
    for (const auto& [name, set] : gMeasured) {
        const std::string label = name;
        ms_hrtf* captured = set;
        run(label.c_str(), [label, captured] { testHeadModelCues(label.c_str(), captured, false); });
    }
    run("diffuse field", testDiffuseField);
    run("symmetry", testSymmetry);
    run("no allocation", testNoAllocationOnRenderSide);
    run("determinism", testDeterminism);
    run("doppler", testDoppler);
    run("distance and air", testDistanceAndAir);
    run("near field", testNearField);
    run("ambisonics", testAmbisonics);
    run("speakers", testSpeakers);
    run("device channels", testDeviceChannels);
    run("crosstalk", testCrosstalk);
    run("moving source", testMovingSourceStaysClean);
    run("voice overflow", testVoiceOverflow);
    run("master gain", testMasterGain);
    run("buses reach everything", testBusesReachEverything);
    run("commit is one batch", testCommitIsOneBatch);
    run("sends hear the source", testSendsHearTheSource);
    run("bad values stay local", testBadValuesStayLocal);
    run("demotion keeps tails", testDemotionKeepsTails);
    run("dry signal ramps", testDrySignalRamps);
    run("low pitch", testLowPitch);
    run("workgroup join", testWorkgroupJoinOnce);
    run("HRTF rate keeps level", testHrtfRateKeepsLevel);
    run("stored delays", testStoredDelaysSurviveImport);
    run("ear order", testEarOrderComesFromTheData);
    run("import rate keeps level", testImportRateKeepsLevel);
    run("scheduling and features", testSchedulingAndFeatures);
    run("block sizes", testBlockSizesAgree);
    run("other rates", testOtherRates);
    run("source shaping", testSourceShaping);
    run("files and SOFA", testFilesAndSofa);
    run("lock-free boundary", testConcurrentCommits);
    run("budget", testBudget);
    for (const auto& [name, set] : gMeasured) ms_hrtf_release(set);
    ms_hrtf_release(gSphere);
    std::printf("%d checks, %d failures\n", gChecks, gFailures);
    return gFailures ? 1 : 0;
}
