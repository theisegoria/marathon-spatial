#include "marathon_spatial/marathon_spatial_acoustics.h"

#include "acoustics/AcousticSimulator.hpp"
#include "acoustics/LateReverb.hpp"
#include "acoustics/MetalTraceBackend.hpp"

#include <algorithm>
#include <cstring>
#include <memory>

using namespace mspatial::acoustics;

struct ms_acoustics {
    Scene scene;
    std::unique_ptr<TraceBackend> backend;
    std::unique_ptr<Simulator> simulator;
    LateReverbNetwork reverb;
    float sampleRate = 48000.f;
    bool gpu = false;
    bool reverbReady = false;
};

namespace {

static_assert(MS_ACOUSTIC_BANDS == kBandCount && MS_ACOUSTIC_MAX_SOURCES == kMaxSources && MS_ACOUSTIC_MAX_ROOMS == kMaxRooms);
static_assert(MS_MAX_REFLECTION_TAPS == kMaxPathsPerSource, "a slot here is a slot in the engine");
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
void ms_acoustics_tick(ms_acoustics* a, float seconds) { a->simulator->tick(seconds); }

uint32_t ms_acoustics_source_taps(const ms_acoustics* a, uint32_t source, ms_reflection_tap* taps) {
    const SourceAcoustics* view = viewOf(a, source);
    uint32_t count = 0;
    for (uint32_t i = 0; i < MS_MAX_REFLECTION_TAPS; ++i) {
        ms_reflection_tap tap{};
        if (view && view->active) {
            const ReflectionPath& path = view->reflections[i];
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
    a->reverb.prepare(config);
    a->reverbReady = true;
    return MS_OK;
}

void ms_acoustics_reverb_process(ms_acoustics* a, const float* const* inputs, float* const* left, float* const* right, uint32_t frames) {
    if (!a->reverbReady) return;
    a->reverb.apply(a->simulator->snapshots().read());
    a->reverb.process(inputs, left, right, frames);
}

}  // extern "C"
