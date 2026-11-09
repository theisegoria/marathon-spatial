/*
 * Marathon Spatial geometric acoustics, C API.
 *
 * Scene as plain data: triangles, materials in octave bands, rooms, portals and moving instances. Hardware ray traced early reflections,
 * occlusion and transmission, room measurement, portal paths and coupled late reverberation. One control thread makes every call except
 * ms_acoustics_reverb_process, which belongs to one render thread, and never waits or allocates.
 *
 * Each frame: move instances and sources, call ms_acoustics_tick, then for each source call ms_acoustics_source_taps and hand the result to
 * ms_source_set_reflections and the source's eq values. Nothing here opens an audio device.
 *
 * The listener's room gets a hybrid late reverb. Paths traced from the listener measure the room's energy field: how much reflected energy
 * arrives, from where and when, per octave band. That field becomes a first-order ambisonic impulse response which convolves the room's send
 * for the first convolution_seconds, and the room's feedback network takes over from there at the level the field reached. The early part
 * stays with the exact image-source taps. Other rooms use their feedback networks alone.
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
    float convolution_seconds;  /* Length of the listener room's convolution reverb. 0 turns it off. Default 1.5. */
    uint32_t energy_paths;      /* Paths tracing the listener's energy field, one segment each per tick. Default 384. 0 leaves only a bake's field. */
    float baked_prior_paths;    /* How many traced paths a bake's field counts as when the live field refines it. Default 4096. */
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

/* The listener's traced energy field and the convolution reverb built from it. */
typedef struct ms_acoustic_energy_view {
    uint32_t active;                   /* A response from the field is in use. */
    uint32_t room;                     /* The room it belongs to, or MS_ACOUSTIC_NONE. */
    float paths;                       /* Effective number of traced paths behind it; old evidence fades over a few seconds. */
    float seconds;                     /* Length of the convolution part. */
    float t60[MS_ACOUSTIC_BANDS];      /* Fitted to the traced field; 0 where it holds too little to fit. */
    float energy[MS_ACOUSTIC_BANDS];   /* Reflected energy of an impulse in the field, relative to the direct sound at 1 m. */
    uint32_t baked;                    /* The field includes a bake's probes. paths then counts them at baked_prior_paths. */
} ms_acoustic_energy_view;

/*
 * Baked propagation for the static part of a scene. A bake places probes through the rooms and stores, for each, the energy field a listener
 * there would measure, and for every pair of probes that see each other, the shortest route between them. Set on an ms_acoustics, it gives the
 * listener's room a converged reverb the moment the listener arrives, which the live field then refines, and routes sound round any number of
 * corners where edge diffraction stops at two. Rays still check every route against the live scene, so a door closed since the bake blocks it.
 */
typedef struct ms_acoustic_bake ms_acoustic_bake;

typedef struct ms_acoustic_bake_config {
    float probe_spacing;         /* Metres between grid probes. Default 2. */
    float clearance;             /* Grid probes closer than this to a surface are dropped. Default 0.3. */
    const ms_vec3* probes;       /* Your own probe positions instead of a grid, or NULL. */
    uint32_t probe_count;
    uint32_t reverb;             /* Bake each probe's energy field. Default 1. */
    uint32_t pathing;            /* Bake routes between probes. Default 1. */
    uint32_t paths_per_probe;    /* Energy paths traced from each probe. Default 4096. */
    uint32_t threads;            /* Default 4. */
} ms_acoustic_bake_config;

void ms_acoustic_bake_config_default(ms_acoustic_bake_config* config);
/* Bakes the scene as it stands, with the field length and air of this ms_acoustics. Slow: seconds to minutes. Call on the control thread;
   progress, when not NULL, is called on it with the share done, 0 to 1. Returns NULL on failure, with the reason in status. */
ms_acoustic_bake* ms_acoustics_bake(const ms_acoustics* acoustics, const ms_acoustic_bake_config* config, void (*progress)(float done, void* user),
                                    void* user, ms_status* status);
ms_status ms_acoustic_bake_save(const ms_acoustic_bake* bake, const char* path);
/* MS_ERR_NOT_FOUND when the file cannot be opened, MS_ERR_BAD_FILE when it is not a whole bake of this version. */
ms_acoustic_bake* ms_acoustic_bake_load(const char* path, ms_status* status);
void ms_acoustic_bake_destroy(ms_acoustic_bake* bake);
uint32_t ms_acoustic_bake_probe_count(const ms_acoustic_bake* bake);
/* Nonzero when the scene's enabled geometry, transforms and materials are what the bake was made from. */
uint32_t ms_acoustic_bake_matches(const ms_acoustic_bake* bake, const ms_acoustics* acoustics);
/* Uses the bake from now on, or none with NULL. It keeps its own reference, so the bake may be destroyed at once. */
void ms_acoustics_set_bake(ms_acoustics* acoustics, const ms_acoustic_bake* bake);

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
/* Which way the listener faces, as for ms_listener_set_pose. Only the stereo decoding of the convolution reverb uses it; the ambisonic
   output of ms_acoustics_reverb_process_ambisonic stays in the world frame. Default: forward -Z, up +Y. */
void ms_acoustics_set_listener_orientation(ms_acoustics* acoustics, ms_vec3 forward, ms_vec3 up);

/* Takes finished rays, refreshes everything analytic, publishes, and queues the next rays. Never waits for the GPU. */
void ms_acoustics_tick(ms_acoustics* acoustics, float seconds);

/* Reads the latest tick. taps must hold MS_MAX_REFLECTION_TAPS; empty slots have zero gain and slots keep their paths, as
   ms_source_set_reflections wants. Returns the count to pass on. */
uint32_t ms_acoustics_source_taps(const ms_acoustics* acoustics, uint32_t source, ms_reflection_tap* taps);
void ms_acoustics_source_view_get(const ms_acoustics* acoustics, uint32_t source, ms_acoustic_source_view* view);
/* paths must hold MS_ACOUSTIC_MAX_PORTAL_PATHS. */
void ms_acoustics_source_portal_paths(const ms_acoustics* acoustics, uint32_t source, ms_acoustic_portal_path* paths);
void ms_acoustics_room_view_get(const ms_acoustics* acoustics, uint32_t room, ms_acoustic_room_view* view);
void ms_acoustics_energy_view_get(const ms_acoustics* acoustics, ms_acoustic_energy_view* view);
/* The room of a source, or the listener with MS_ACOUSTIC_NONE as source. */
uint32_t ms_acoustics_room_of(const ms_acoustics* acoustics, uint32_t source);

/* Late reverb: one network per room. Allocates; call on the control thread once rooms exist and with no render running.
   Tail delays follow each room's mean free path at this moment and stay fixed until the next prepare. */
ms_status ms_acoustics_reverb_prepare(ms_acoustics* acoustics);
/* Render thread. inputs[room] is the mono send into that room or NULL; left[room] and right[room] are outputs or NULL. Arrays hold one entry per
   room. Output is pressure relative to free field at 1 m, so the send level is 1 for a physical result. */
void ms_acoustics_reverb_process(ms_acoustics* acoustics, const float* const* inputs, float* const* left, float* const* right, uint32_t frames);
/* As ms_acoustics_reverb_process, but the listener room's convolution reverb leaves as first-order ambisonics in ambisonic[0..3] (ACN/SN3D
   AmbiX in the world frame, for a world-locked ms_bed, which decodes it through the HRTF and follows the head) instead of being decoded into
   that room's left and right. ambisonic entries may be NULL. */
void ms_acoustics_reverb_process_ambisonic(ms_acoustics* acoustics, const float* const* inputs, float* const* left, float* const* right,
                                           float* const* ambisonic, uint32_t frames);

/* Binaural late reverberation. Each room network's left and right become mixes of its two orthogonal outputs whose interaural coherence
   follows coherence[band] in six bands centred on 141, 283, 566, 1131, 2263 and 4525 Hz (crossovers 200, 400, 800, 1600, 3200 Hz), at
   unchanged level. Give it the listener's HRTF set's diffuse-field coherence (ms_hrtf_diffuse_coherence at those centres), so the tail
   reaches the ears as a diffuse field would. NULL returns to two uncorrelated outputs. Any thread; heard from the next render chunk. */
#define MS_LATE_COHERENCE_BANDS 6
extern const float ms_late_coherence_hz[MS_LATE_COHERENCE_BANDS];
void ms_acoustics_set_late_coherence(ms_acoustics* acoustics, const float* coherence);

#ifdef __cplusplus
}
#endif

#endif
