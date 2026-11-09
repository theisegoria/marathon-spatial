#include "marathon_spatial/marathon_spatial_acoustics.h"

#include "acoustics/AcousticSimulator.hpp"
#include "acoustics/Bake.hpp"
#include "acoustics/Convolver.hpp"
#include "acoustics/LateReverb.hpp"
#include "acoustics/MetalTraceBackend.hpp"
#include "acoustics/ReverbIr.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <functional>
#include <memory>
#include <numbers>

using namespace mspatial::acoustics;

namespace {

constexpr std::uint32_t kConvolutionBlock = 256;
constexpr std::uint32_t kRenderChunk = 512;
constexpr float kIrInterval = 0.5f;     // Seconds between responses built from the field.
constexpr double kMinimumPaths = 64.0;  // Evidence below this is not worth a response.
// Virtual microphones for stereo: a + (1 - a) cos(theta) with a = 1 / (1 + sqrt 3) leaves a diffuse field uncorrelated between the ears.
constexpr float kMicOmni = 0.3660254f;

// What the render thread needs to run the hybrid, published by the control thread with each new response.
struct HybridParams {
    std::uint32_t room = kNone;      // Convolution input and output room, kNone when off.
    std::uint32_t delay = 0;         // The network's delayed send, in samples.
    Bands gainDb{};                  // And its level, per band.
    Vec3 forward = vec3(0, 0, -1);
    Vec3 up = vec3(0, 1, 0);
};

}  // namespace

struct ms_acoustic_bake {
    std::shared_ptr<const AcousticBake> bake;
};

struct ms_acoustics {
    Scene scene;
    std::unique_ptr<TraceBackend> backend;
    std::unique_ptr<Simulator> simulator;
    LateReverbNetwork reverb;
    float sampleRate = 48000.f;
    float convolutionSeconds = 0.f;
    bool gpu = false;
    bool reverbReady = false;

    // Control thread.
    Bands air = airAttenuation();
    std::shared_ptr<const AcousticBake> bake;
    double bakedPrior = 4096.0;
    bool bakedInField = false;
    ReverbIrSynth synth;
    EnergyField field, live, baked;
    std::vector<float> irStore;
    float sinceResponse = 0.f;
    bool responseActive = false;
    std::uint32_t responseRoom = kNone;
    HybridParams params;

    // Render thread.
    PartitionedConvolver convolver;
    SnapshotExchange<HybridParams> hybrid;
    std::uint32_t convolverRoom = kNone;
    std::vector<float> renderScratch;  // Four ambisonic channels of kRenderChunk.
};

namespace {

static_assert(MS_ACOUSTIC_BANDS == kBandCount && MS_ACOUSTIC_MAX_SOURCES == kMaxSources && MS_ACOUSTIC_MAX_ROOMS == kMaxRooms);
static_assert(MS_MAX_REFLECTION_TAPS == kMaxPathsPerSource + kMaxDiffractionPaths, "a slot here is a slot in the engine");
static_assert(MS_ACOUSTIC_MAX_PORTAL_PATHS == kMaxPortalPaths);

Vec3 toVec3(ms_vec3 v) { return vec3(v.x, v.y, v.z); }
ms_vec3 fromVec3(Vec3 v) { return {v.x, v.y, v.z}; }

Bands toBands(const float (&values)[MS_ACOUSTIC_BANDS]) {
    Bands out;
    std::copy_n(values, kBandCount, out.begin());
    return out;
}

Material toMaterial(const ms_acoustic_material& m) { return {toBands(m.absorption), toBands(m.scattering), toBands(m.transmission)}; }

Transform toTransform(const ms_acoustic_transform* t) {
    Transform out;
    if (t) std::memcpy(out.rows, t->rows, sizeof out.rows);
    return out;
}

PortalDesc toPortal(const ms_acoustic_portal& p) {
    PortalDesc out;
    out.roomA = p.room_a;
    out.roomB = p.room_b;
    out.center = toVec3(p.center);
    out.normal = toVec3(p.normal);
    out.halfU = toVec3(p.half_u);
    out.halfV = toVec3(p.half_v);
    out.halfThickness = p.half_thickness;
    return out;
}

// The engine's three bands: low below 400 Hz, mid to 2.5 kHz, high above.
void threeBands(const Bands& in, float out[3]) {
    out[0] = (in[0] + in[1] + in[2]) / 3.f;
    out[1] = (in[3] + in[4] + in[5]) / 3.f;
    out[2] = (in[6] + in[7]) / 2.f;
}

const SourceAcoustics* viewOf(const ms_acoustics* a, std::uint32_t source) {
    return source < kMaxSources && a->simulator->lastPublished() ? &a->simulator->lastPublished()->sources[source] : nullptr;
}

}  // namespace

extern "C" {

void ms_acoustics_config_default(ms_acoustics_config* config) {
    *config = {};
    config->sample_rate = 48000.f;
    config->temperature_celsius = 20.f;
    config->humidity_percent = 50.f;
    config->convolution_seconds = 1.5f;
    config->energy_paths = 384;
    config->baked_prior_paths = 4096.f;
}

ms_acoustics* ms_acoustics_create(const ms_acoustics_config* config) {
    ms_acoustics_config defaults;
    ms_acoustics_config_default(&defaults);
    if (!config) config = &defaults;
    auto a = std::make_unique<ms_acoustics>();
    if (config->sample_rate > 0.f) a->sampleRate = config->sample_rate;
    if (config->metallib_path) {
        MetalTraceConfig metal;
        metal.libraryPath = config->metallib_path;
        std::string error;
        a->backend = makeMetalTraceBackend(metal, &error);
        a->gpu = a->backend != nullptr;
    }
    if (!a->backend) a->backend = makeCpuTraceBackend();
    SimulatorConfig simulation;
    simulation.air = airAttenuation(config->temperature_celsius, config->humidity_percent);
    simulation.formula = config->use_sabine ? DecayFormula::Sabine : DecayFormula::Eyring;
    a->air = simulation.air;
    a->bakedPrior = std::max(0.f, config->baked_prior_paths);
    a->convolutionSeconds = std::clamp(config->convolution_seconds, 0.f, 4.f);
    // With no live paths the convolution still runs from a bake's field, once one is set.
    simulation.energy.paths = a->convolutionSeconds > 0.f ? config->energy_paths : 0;
    simulation.energy.seconds = std::max(a->convolutionSeconds, 0.01f);
    a->simulator = std::make_unique<Simulator>(a->scene, *a->backend, simulation);
    return a.release();
}

void ms_acoustics_destroy(ms_acoustics* acoustics) { delete acoustics; }
uint32_t ms_acoustics_uses_gpu(const ms_acoustics* acoustics) { return acoustics->gpu; }

uint32_t ms_acoustics_material_add(ms_acoustics* a, const ms_acoustic_material* material) { return a->scene.addMaterial(toMaterial(*material)); }
void ms_acoustics_material_set(ms_acoustics* a, uint32_t material, const ms_acoustic_material* value) { a->scene.setMaterial(material, toMaterial(*value)); }

uint32_t ms_acoustics_mesh_add(ms_acoustics* a, const float* positions, uint32_t vertexCount, const uint32_t* indices, uint32_t indexCount,
                               const uint32_t* triangleMaterials) {
    return a->scene.addMesh({{positions, std::size_t(vertexCount) * 3}, {indices, indexCount}, {triangleMaterials, indexCount / 3}});
}
void ms_acoustics_mesh_remove(ms_acoustics* a, uint32_t mesh) { a->scene.removeMesh(mesh); }

uint32_t ms_acoustics_instance_add(ms_acoustics* a, uint32_t mesh, const ms_acoustic_transform* transform) {
    return a->scene.addInstance(mesh, toTransform(transform));
}
void ms_acoustics_instance_set_transform(ms_acoustics* a, uint32_t instance, const ms_acoustic_transform* t) { a->scene.setInstanceTransform(instance, toTransform(t)); }
void ms_acoustics_instance_set_enabled(ms_acoustics* a, uint32_t instance, int enabled) { a->scene.setInstanceEnabled(instance, enabled != 0); }
void ms_acoustics_instance_remove(ms_acoustics* a, uint32_t instance) { a->scene.removeInstance(instance); }

uint32_t ms_acoustics_room_add(ms_acoustics* a, const ms_acoustic_room* room) {
    RoomDesc desc;
    for (uint32_t i = 0; i < room->bound_count; ++i) desc.bounds.push_back({toVec3(room->bounds[i].min), toVec3(room->bounds[i].max)});
    desc.probe = toVec3(room->probe);
    desc.volume = room->volume;
    desc.surfaceArea = room->surface_area;
    desc.absorptionArea = toBands(room->absorption_area);
    return a->scene.addRoom(desc);
}
uint32_t ms_acoustics_portal_add(ms_acoustics* a, const ms_acoustic_portal* portal) { return a->scene.addPortal(toPortal(*portal)); }
void ms_acoustics_portal_set(ms_acoustics* a, uint32_t portal, const ms_acoustic_portal* value) { a->scene.setPortal(portal, toPortal(*value)); }

uint32_t ms_acoustics_source_add(ms_acoustics* a, ms_vec3 position, float radius) { return a->simulator->addSource({toVec3(position), radius}); }
void ms_acoustics_source_set(ms_acoustics* a, uint32_t source, ms_vec3 position, float radius) { a->simulator->setSource(source, {toVec3(position), radius}); }
void ms_acoustics_source_remove(ms_acoustics* a, uint32_t source) { a->simulator->removeSource(source); }
void ms_acoustics_set_listener(ms_acoustics* a, ms_vec3 position) { a->simulator->setListener(toVec3(position)); }

void ms_acoustics_set_listener_orientation(ms_acoustics* a, ms_vec3 forward, ms_vec3 up) {
    const Vec3 f = toVec3(forward), u = toVec3(up);
    if (length(f) < 1e-6f || length(cross(f, u)) < 1e-6f) return;
    a->params.forward = normalize(f);
    a->params.up = normalize(u);
    a->hybrid.writeSlot() = a->params;
    a->hybrid.publish();
}

namespace {

// A new response from the listener's field, or silence when there is not enough of it. Control thread.
void buildResponse(ms_acoustics* a) {
    const AcousticSnapshot* snapshot = a->simulator->lastPublished();
    const EnergyTracer& tracer = a->simulator->energyTracer();
    tracer.read(a->live);
    // A bake's probes in the listener's room give the field a converged start, worth bakedPrior paths, which the live paths then refine.
    const float binSeconds = tracer.config().binSeconds;
    const auto bins = std::uint32_t(std::ceil(a->convolutionSeconds / binSeconds));
    const std::uint32_t listenerRoom = snapshot ? snapshot->listenerRoom : kNone;
    a->bakedInField = a->bake && a->bake->hasReverb() && listenerRoom != kNone && a->bakedPrior > 0.0 &&
                      a->bake->fieldAt(a->simulator->listener(), listenerRoom, bins, binSeconds, a->baked);
    if (a->bakedInField) {
        if (a->live.room == listenerRoom && a->live.bins == bins && a->live.paths > 0.0) {
            blendFields(a->live, a->baked, a->bakedPrior, a->field);
        } else {
            a->field = a->baked;
            a->field.paths = a->bakedPrior;
        }
    } else {
        a->field = a->live;
    }
    const std::uint32_t room = a->field.room;
    const bool usable = snapshot && room != kNone && room < snapshot->roomCount && snapshot->rooms[room].active && a->field.paths >= kMinimumPaths;
    if (!usable) {
        if (a->responseActive && a->convolver.stageSilence()) {
            a->responseActive = false;
            a->params.room = kNone;
            a->hybrid.writeSlot() = a->params;
            a->hybrid.publish();
        }
        return;
    }
    const RoomAcoustics& acoustics = snapshot->rooms[room];
    const float rate = a->sampleRate, length = a->convolutionSeconds;
    // The image-source taps own the early part; the response rises over about the mixing time, sqrt(V) milliseconds.
    const float mixing = std::clamp(0.001f * std::sqrt(std::max(acoustics.volume, 1.f)), 0.02f, 0.12f);
    ReverbIrShape shape;
    shape.advance = a->convolver.latency();
    shape.fadeInStart = std::max(1.2f * float(shape.advance) / rate, 0.5f * mixing);
    shape.fadeInEnd = std::max(shape.fadeInStart + 0.005f, 1.5f * mixing);
    // The network builds up over its first few line delays, about a mean free path each, so the response fades out over that much: the two
    // cross over in energy with neither a gap nor a doubling.
    const float handOver = std::clamp(3.f * acoustics.meanFreePath / kDefaultSpeedOfSound, 0.02f, 0.15f);
    shape.fadeOutStart = std::max(shape.fadeInEnd, length - handOver);
    shape.fadeOutEnd = length;
    float* channels[kAmbisonicChannels];
    for (std::uint32_t c = 0; c < kAmbisonicChannels; ++c) channels[c] = a->irStore.data() + std::size_t(c) * a->synth.maxSamples();
    const std::uint32_t samples = a->synth.build(a->field, shape, channels);
    if (!a->convolver.stage(channels, samples)) return;  // The previous one is still being picked up: try again next tick.
    a->responseActive = true;

    // The network starts at the start of the fade, at the level the field holds there. Its own initial energy rate is 4 pi c / V per
    // second relative to the direct sound at 1 m (16 pi / A of energy decaying at c A / 4 V).
    const EnergyField& smoothed = a->synth.smoothed();
    const auto bin = std::min<std::uint32_t>(std::uint32_t(shape.fadeOutStart / smoothed.binSeconds), smoothed.bins - 1);
    const float networkRate = 4.f * std::numbers::pi_v<float> * kDefaultSpeedOfSound / std::max(acoustics.volume, 1.f);
    for (std::uint32_t b = 0; b < kBandCount; ++b) {
        const float traced = smoothed.band(b, bin) / smoothed.binSeconds;
        const float natural = -4.342945f * acoustics.decayRate[b] * shape.fadeOutStart;
        a->params.gainDb[b] = std::clamp(traced > 0.f ? 10.f * std::log10(traced / networkRate) : natural, -150.f, 12.f);
    }
    a->params.room = room;
    a->params.delay = std::uint32_t(shape.fadeOutStart * rate);
    a->hybrid.writeSlot() = a->params;
    a->hybrid.publish();
}

}  // namespace

void ms_acoustics_tick(ms_acoustics* a, float seconds) {
    a->simulator->tick(seconds);
    if (!a->reverbReady || a->convolutionSeconds <= 0.f) return;
    a->sinceResponse += seconds;
    if (a->sinceResponse < kIrInterval) return;
    a->sinceResponse = 0.f;
    buildResponse(a);
}

void ms_acoustics_energy_view_get(const ms_acoustics* a, ms_acoustic_energy_view* out) {
    *out = {};
    out->room = MS_ACOUSTIC_NONE;
    EnergyField field;
    if (a->bakedInField)
        field = a->field;
    else
        a->simulator->energyTracer().read(field);
    out->baked = a->bakedInField;
    out->active = a->responseActive;
    out->room = a->responseActive ? a->params.room : field.room;
    out->paths = float(field.paths);
    out->seconds = a->convolutionSeconds;
    const Bands t60 = field.t60();
    std::copy(t60.begin(), t60.end(), out->t60);
    for (std::uint32_t b = 0; b < kBandCount; ++b) {
        float sum = 0.f;
        for (std::uint32_t i = 0; i < field.bins; ++i) sum += field.band(b, i);
        out->energy[b] = sum;
    }
}

uint32_t ms_acoustics_source_taps(const ms_acoustics* a, uint32_t source, ms_reflection_tap* taps) {
    const SourceAcoustics* view = viewOf(a, source);
    uint32_t count = 0;
    for (uint32_t i = 0; i < MS_MAX_REFLECTION_TAPS; ++i) {
        ms_reflection_tap tap{};
        if (view && view->active) {
            // Early reflections first, then the diffracted paths in the last slots.
            const ReflectionPath& path = i < kMaxPathsPerSource ? view->reflections[i] : view->diffraction[i - kMaxPathsPerSource];
            if (path.id != 0 && path.excessSeconds <= 0.17f) {
                tap.delay_seconds = path.excessSeconds;
                tap.direction = fromVec3(path.arrival);
                threeBands(path.gain, tap.gain);
                count = i + 1;
            }
        }
        taps[i] = tap;
    }
    return count;
}

void ms_acoustics_source_view_get(const ms_acoustics* a, uint32_t source, ms_acoustic_source_view* out) {
    *out = {};
    out->room = MS_ACOUSTIC_NONE;
    const SourceAcoustics* view = viewOf(a, source);
    if (!view || !view->active) return;
    out->active = 1;
    out->room = view->room;
    out->visibility = view->visibility;
    float eq[3];
    threeBands(view->directGain, eq);
    const auto db = [](float gain) { return 20.f * std::log10(std::max(gain, 1e-4f)); };
    out->eq_low_db = db(eq[0]);
    out->eq_mid_db = db(eq[1]);
    out->eq_high_db = db(eq[2]);
}

void ms_acoustics_source_portal_paths(const ms_acoustics* a, uint32_t source, ms_acoustic_portal_path* paths) {
    const SourceAcoustics* view = viewOf(a, source);
    for (uint32_t i = 0; i < MS_ACOUSTIC_MAX_PORTAL_PATHS; ++i) {
        paths[i] = {};
        if (!view || !view->active || view->portalPaths[i].portals == 0) continue;
        const PortalPath& path = view->portalPaths[i];
        threeBands(path.gain, paths[i].gain);
        paths[i].delay_seconds = path.excessSeconds;
        paths[i].direction = fromVec3(path.arrival);
        paths[i].portals = path.portals;
    }
}

void ms_acoustics_room_view_get(const ms_acoustics* a, uint32_t room, ms_acoustic_room_view* out) {
    *out = {};
    const AcousticSnapshot* snapshot = a->simulator->lastPublished();
    if (!snapshot || room >= snapshot->roomCount || !snapshot->rooms[room].active) return;
    const RoomAcoustics& view = snapshot->rooms[room];
    out->active = 1;
    out->volume = view.volume;
    out->surface_area = view.surfaceArea;
    std::copy(view.t60.begin(), view.t60.end(), out->t60);
}

uint32_t ms_acoustics_room_of(const ms_acoustics* a, uint32_t source) {
    const AcousticSnapshot* snapshot = a->simulator->lastPublished();
    if (!snapshot) return MS_ACOUSTIC_NONE;
    if (source == MS_ACOUSTIC_NONE) return snapshot->listenerRoom;
    const SourceAcoustics* view = viewOf(a, source);
    return view && view->active ? view->room : MS_ACOUSTIC_NONE;
}

ms_status ms_acoustics_reverb_prepare(ms_acoustics* a) {
    LateReverbConfig config;
    config.sampleRate = a->sampleRate;
    config.rooms = std::min<std::uint32_t>(std::uint32_t(a->scene.rooms().size()), kMaxRooms);
    config.links = std::min<std::uint32_t>(std::uint32_t(a->scene.portals().size()) * 2, kMaxRoomLinks);
    if (const AcousticSnapshot* snapshot = a->simulator->lastPublished()) {
        for (std::uint32_t r = 0; r < config.rooms; ++r) config.meanDelaySeconds.push_back(snapshot->rooms[r].meanFreePath / kDefaultSpeedOfSound);
    }
    if (a->convolutionSeconds > 0.f) {
        const auto samples = std::uint32_t(a->convolutionSeconds * a->sampleRate);
        config.maxHybridDelay = samples;
        a->convolver.prepare(kConvolutionBlock, samples, kAmbisonicChannels);
        a->synth.prepare(a->sampleRate, samples);
        a->irStore.assign(std::size_t(kAmbisonicChannels) * a->synth.maxSamples(), 0.f);
        a->renderScratch.assign(std::size_t(kAmbisonicChannels) * kRenderChunk, 0.f);
    }
    a->reverb.prepare(config);
    a->responseActive = false;
    a->convolverRoom = kNone;
    a->params.room = kNone;
    a->hybrid.writeSlot() = a->params;
    a->hybrid.publish();
    a->reverbReady = true;
    return MS_OK;
}

namespace {

void processHybrid(ms_acoustics* a, const float* const* inputs, float* const* left, float* const* right, float* const* ambisonic, uint32_t frames) {
    if (!a->reverbReady) return;
    a->reverb.apply(a->simulator->snapshots().read());
    const HybridParams& params = a->hybrid.read();
    if (!a->convolver.prepared()) {
        a->reverb.process(inputs, left, right, frames);
        return;
    }
    // The convolver keeps the room it last had until a response for another room arrives, so a fading tail stays where it was.
    if (params.room != kNone && params.room != a->convolverRoom) {
        if (a->convolverRoom != kNone) a->convolver.reset();
        a->convolverRoom = params.room;
    }
    a->reverb.setHybrid(params.room, params.delay, params.gainDb);
    const std::uint32_t rooms = std::min<std::uint32_t>(std::uint32_t(a->scene.rooms().size()), kMaxRooms);
    const float* inputAt[kMaxRooms];
    float* leftAt[kMaxRooms];
    float* rightAt[kMaxRooms];
    float* foa[kAmbisonicChannels];
    for (std::uint32_t c = 0; c < kAmbisonicChannels; ++c) foa[c] = a->renderScratch.data() + std::size_t(c) * kRenderChunk;
    const Vec3 right3 = normalize(cross(params.forward, params.up));
    const float micGain = 1.f / std::sqrt(kMicOmni * kMicOmni + (1.f - kMicOmni) * (1.f - kMicOmni) / 3.f);
    // A virtual microphone pointing along world direction v sees AmbiX components (X, Y, Z) = (-v.z, -v.x, v.y).
    const auto weights = [&](Vec3 v, float (&w)[kAmbisonicChannels]) {
        w[0] = kMicOmni * micGain;
        w[1] = (1.f - kMicOmni) * micGain * -v.x;
        w[2] = (1.f - kMicOmni) * micGain * v.y;
        w[3] = (1.f - kMicOmni) * micGain * -v.z;
    };
    float leftWeights[kAmbisonicChannels], rightWeights[kAmbisonicChannels];
    weights(-right3, leftWeights);
    weights(right3, rightWeights);
    std::uint32_t done = 0;
    while (done < frames) {
        const std::uint32_t chunk = std::min(frames - done, kRenderChunk);
        for (std::uint32_t r = 0; r < rooms; ++r) {
            inputAt[r] = inputs && inputs[r] ? inputs[r] + done : nullptr;
            leftAt[r] = left && left[r] ? left[r] + done : nullptr;
            rightAt[r] = right && right[r] ? right[r] + done : nullptr;
        }
        a->reverb.process(inputAt, leftAt, rightAt, chunk);
        const std::uint32_t room = a->convolverRoom;
        if (room == kNone || room >= rooms) {
            if (ambisonic)
                for (std::uint32_t c = 0; c < kAmbisonicChannels; ++c)
                    if (ambisonic[c]) std::fill_n(ambisonic[c] + done, chunk, 0.f);
            done += chunk;
            continue;
        }
        a->convolver.process(inputAt[room], foa, chunk);
        if (ambisonic) {
            for (std::uint32_t c = 0; c < kAmbisonicChannels; ++c)
                if (ambisonic[c]) std::copy_n(foa[c], chunk, ambisonic[c] + done);
        } else {
            for (std::uint32_t c = 0; c < kAmbisonicChannels; ++c) {
                if (leftAt[room]) vDSP_vsma(foa[c], 1, &leftWeights[c], leftAt[room], 1, leftAt[room], 1, chunk);
                if (rightAt[room]) vDSP_vsma(foa[c], 1, &rightWeights[c], rightAt[room], 1, rightAt[room], 1, chunk);
            }
        }
        done += chunk;
    }
}

}  // namespace

const float ms_late_coherence_hz[MS_LATE_COHERENCE_BANDS] = {141.4f, 282.8f, 565.7f, 1131.4f, 2262.7f, 4525.5f};

void ms_acoustics_set_late_coherence(ms_acoustics* a, const float* coherence) {
    if (a) a->reverb.setCoherence(coherence);
}

void ms_acoustics_reverb_process(ms_acoustics* a, const float* const* inputs, float* const* left, float* const* right, uint32_t frames) {
    processHybrid(a, inputs, left, right, nullptr, frames);
}

void ms_acoustics_reverb_process_ambisonic(ms_acoustics* a, const float* const* inputs, float* const* left, float* const* right,
                                           float* const* ambisonic, uint32_t frames) {
    static float* const none[kAmbisonicChannels] = {};
    processHybrid(a, inputs, left, right, ambisonic ? ambisonic : none, frames);
}

void ms_acoustic_bake_config_default(ms_acoustic_bake_config* config) {
    *config = {};
    config->probe_spacing = 2.f;
    config->clearance = 0.3f;
    config->reverb = 1;
    config->pathing = 1;
    config->paths_per_probe = 4096;
    config->threads = 4;
}

ms_acoustic_bake* ms_acoustics_bake(const ms_acoustics* a, const ms_acoustic_bake_config* config, void (*progress)(float, void*), void* user,
                                    ms_status* status) {
    ms_acoustic_bake_config defaults;
    ms_acoustic_bake_config_default(&defaults);
    if (!config) config = &defaults;
    const auto report = [&](ms_status s) {
        if (status) *status = s;
    };
    if (!a || (config->probe_count > 0 && !config->probes)) {
        report(MS_ERR_INVALID_ARGUMENT);
        return nullptr;
    }
    BakeConfig bake;
    bake.spacing = config->probe_spacing;
    bake.clearance = config->clearance;
    for (uint32_t i = 0; i < config->probe_count; ++i) bake.probes.push_back(toVec3(config->probes[i]));
    bake.reverb = config->reverb != 0;
    bake.pathing = config->pathing != 0;
    bake.pathsPerProbe = config->paths_per_probe;
    bake.threads = config->threads;
    bake.energy.seconds = a->convolutionSeconds > 0.f ? a->convolutionSeconds : 1.5f;
    bake.energy.air = a->air;
    auto out = std::make_unique<ms_acoustic_bake>();
    auto result = std::make_shared<AcousticBake>();
    std::function<void(float)> callback;
    if (progress) callback = [&](float done) { progress(done, user); };
    if (!bakeScene(a->scene, bake, *result, nullptr, callback)) {
        report(bake.probes.size() > kMaxBakeProbes ? MS_ERR_LIMIT_REACHED : MS_ERR_INVALID_ARGUMENT);
        return nullptr;
    }
    out->bake = std::move(result);
    report(MS_OK);
    return out.release();
}

ms_status ms_acoustic_bake_save(const ms_acoustic_bake* bake, const char* path) {
    if (!bake || !path) return MS_ERR_INVALID_ARGUMENT;
    return bake->bake->save(path) ? MS_OK : MS_ERR_NOT_FOUND;
}

ms_acoustic_bake* ms_acoustic_bake_load(const char* path, ms_status* status) {
    const auto report = [&](ms_status s) {
        if (status) *status = s;
    };
    if (!path) {
        report(MS_ERR_INVALID_ARGUMENT);
        return nullptr;
    }
    if (!std::ifstream(path, std::ios::binary)) {
        report(MS_ERR_NOT_FOUND);
        return nullptr;
    }
    auto result = std::make_shared<AcousticBake>();
    if (!AcousticBake::load(path, *result)) {
        report(MS_ERR_BAD_FILE);
        return nullptr;
    }
    auto out = std::make_unique<ms_acoustic_bake>();
    out->bake = std::move(result);
    report(MS_OK);
    return out.release();
}

void ms_acoustic_bake_destroy(ms_acoustic_bake* bake) { delete bake; }
uint32_t ms_acoustic_bake_probe_count(const ms_acoustic_bake* bake) { return bake ? bake->bake->probeCount() : 0; }
uint32_t ms_acoustic_bake_matches(const ms_acoustic_bake* bake, const ms_acoustics* a) {
    return bake && a && bake->bake->sceneHash == staticSceneHash(a->scene);
}

void ms_acoustics_set_bake(ms_acoustics* a, const ms_acoustic_bake* bake) {
    a->bake = bake ? bake->bake : nullptr;
    a->simulator->setBake(a->bake);
}

}  // extern "C"
