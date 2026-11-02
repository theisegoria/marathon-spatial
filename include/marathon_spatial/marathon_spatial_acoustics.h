/*
 * Marathon Spatial geometric acoustics, C API.
 *
 * Scene as plain data: triangles, materials in octave bands, rooms, portals and moving instances. Hardware ray traced early reflections,
 * occlusion and transmission, room measurement, portal paths and coupled late reverberation. One control thread makes every call except
 * ms_acoustics_reverb_process, which belongs to one render thread, and never waits or allocates.
 *
 * Each frame: move instances and sources, call ms_acoustics_tick, then for each source call ms_acoustics_source_taps and hand the result to
 * ms_source_set_reflections and the source's eq values. Nothing here opens an audio device.
 */
#ifndef MARATHON_SPATIAL_ACOUSTICS_H
#define MARATHON_SPATIAL_ACOUSTICS_H

#include "marathon_spatial.h"
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define MS_ACOUSTIC_BANDS 8u
#define MS_ACOUSTIC_NONE 0xFFFFFFFFu
#define MS_ACOUSTIC_MAX_SOURCES 64u
#define MS_ACOUSTIC_MAX_ROOMS 32u
#define MS_ACOUSTIC_MAX_PORTAL_PATHS 3u

typedef struct ms_acoustics ms_acoustics;

typedef struct ms_acoustics_config {
    const char* metallib_path;  /* acoustic_trace.metallib. NULL or an unusable machine runs the same queries on the calling thread. */
    float sample_rate;          /* For the late reverb. 0 means 48000. */
    float temperature_celsius;  /* Air absorption, ISO 9613-1. */
    float humidity_percent;
    uint32_t use_sabine;        /* Nonzero uses Sabine decay, zero Eyring. */
} ms_acoustics_config;

/* Energy fractions per octave band, 62.5 Hz to 8 kHz. */
typedef struct ms_acoustic_material {
    float absorption[MS_ACOUSTIC_BANDS];
    float scattering[MS_ACOUSTIC_BANDS];
    float transmission[MS_ACOUSTIC_BANDS];
} ms_acoustic_material;

/* Row-major 3 x 4 object-to-world transform. */
typedef struct ms_acoustic_transform {
    float rows[3][4];
} ms_acoustic_transform;

typedef struct ms_acoustic_box {
    ms_vec3 min, max;
} ms_acoustic_box;

/* Leave volume at 0 to measure the room with rays from probe, which follows moving doors and destruction. */
typedef struct ms_acoustic_room {
    const ms_acoustic_box* bounds;
    uint32_t bound_count;
    ms_vec3 probe;
    float volume;
    float surface_area;
    float absorption_area[MS_ACOUSTIC_BANDS];
} ms_acoustic_room;

/* normal points from room_a to room_b; half_u and half_v span the opening. */
typedef struct ms_acoustic_portal {
    uint32_t room_a, room_b;
    ms_vec3 center, normal, half_u, half_v;
    float half_thickness;
} ms_acoustic_portal;

typedef struct ms_acoustic_source_view {
    uint32_t active;
    uint32_t room;
    float visibility;                    /* Fraction of the source disc with a clear line to the listener. */
    float eq_low_db, eq_mid_db, eq_high_db;  /* Occlusion and transmission on the direct path, for ms_source_params. */
} ms_acoustic_source_view;

typedef struct ms_acoustic_portal_path {
    float gain[3];        /* Linear, low, mid, high, relative to free field at 1 m. */
    float delay_seconds;  /* (route - straight line) / speed of sound, as ms_reflection_tap.delay_seconds. */
    ms_vec3 direction;    /* Where the sound seems to come from, seen from the listener. */
    uint32_t portals;     /* 0 for an empty entry. */
} ms_acoustic_portal_path;

typedef struct ms_acoustic_room_view {
    uint32_t active;
    float volume, surface_area;
    float t60[MS_ACOUSTIC_BANDS];
} ms_acoustic_room_view;

void ms_acoustics_config_default(ms_acoustics_config* config);
ms_acoustics* ms_acoustics_create(const ms_acoustics_config* config);
void ms_acoustics_destroy(ms_acoustics* acoustics);
/* Nonzero when rays run on the GPU. */
uint32_t ms_acoustics_uses_gpu(const ms_acoustics* acoustics);

uint32_t ms_acoustics_material_add(ms_acoustics* acoustics, const ms_acoustic_material* material);
void ms_acoustics_material_set(ms_acoustics* acoustics, uint32_t material, const ms_acoustic_material* value);

/* positions holds x,y,z triples; triangle_materials has one id per triangle. */
uint32_t ms_acoustics_mesh_add(ms_acoustics* acoustics, const float* positions, uint32_t vertex_count, const uint32_t* indices, uint32_t index_count,
                               const uint32_t* triangle_materials);
void ms_acoustics_mesh_remove(ms_acoustics* acoustics, uint32_t mesh);

uint32_t ms_acoustics_instance_add(ms_acoustics* acoustics, uint32_t mesh, const ms_acoustic_transform* transform);
void ms_acoustics_instance_set_transform(ms_acoustics* acoustics, uint32_t instance, const ms_acoustic_transform* transform);
void ms_acoustics_instance_set_enabled(ms_acoustics* acoustics, uint32_t instance, int enabled);
void ms_acoustics_instance_remove(ms_acoustics* acoustics, uint32_t instance);

/* Add rooms and portals before the reverb is prepared with ms_acoustics_reverb_prepare. */
uint32_t ms_acoustics_room_add(ms_acoustics* acoustics, const ms_acoustic_room* room);
uint32_t ms_acoustics_portal_add(ms_acoustics* acoustics, const ms_acoustic_portal* portal);
void ms_acoustics_portal_set(ms_acoustics* acoustics, uint32_t portal, const ms_acoustic_portal* value);

uint32_t ms_acoustics_source_add(ms_acoustics* acoustics, ms_vec3 position, float radius);
void ms_acoustics_source_set(ms_acoustics* acoustics, uint32_t source, ms_vec3 position, float radius);
void ms_acoustics_source_remove(ms_acoustics* acoustics, uint32_t source);
void ms_acoustics_set_listener(ms_acoustics* acoustics, ms_vec3 position);

/* Takes finished rays, refreshes everything analytic, publishes, and queues the next rays. Never waits for the GPU. */
void ms_acoustics_tick(ms_acoustics* acoustics, float seconds);

/* Reads the latest tick. taps must hold MS_MAX_REFLECTION_TAPS; empty slots have zero gain and slots keep their paths, as
   ms_source_set_reflections wants. Returns the count to pass on. */
uint32_t ms_acoustics_source_taps(const ms_acoustics* acoustics, uint32_t source, ms_reflection_tap* taps);
void ms_acoustics_source_view_get(const ms_acoustics* acoustics, uint32_t source, ms_acoustic_source_view* view);
/* paths must hold MS_ACOUSTIC_MAX_PORTAL_PATHS. */
void ms_acoustics_source_portal_paths(const ms_acoustics* acoustics, uint32_t source, ms_acoustic_portal_path* paths);
void ms_acoustics_room_view_get(const ms_acoustics* acoustics, uint32_t room, ms_acoustic_room_view* view);
/* The room of a source, or the listener with MS_ACOUSTIC_NONE as source. */
uint32_t ms_acoustics_room_of(const ms_acoustics* acoustics, uint32_t source);

/* Late reverb: one network per room. Allocates; call on the control thread once rooms exist and with no render running.
   Tail delays follow each room's mean free path at this moment and stay fixed until the next prepare. */
ms_status ms_acoustics_reverb_prepare(ms_acoustics* acoustics);
/* Render thread. inputs[room] is the mono send into that room or NULL; left[room] and right[room] are outputs or NULL. Arrays hold one entry per
   room. Output is pressure relative to free field at 1 m, so the send level is 1 for a physical result. */
void ms_acoustics_reverb_process(ms_acoustics* acoustics, const float* const* inputs, float* const* left, float* const* right, uint32_t frames);

#ifdef __cplusplus
}
#endif

#endif
