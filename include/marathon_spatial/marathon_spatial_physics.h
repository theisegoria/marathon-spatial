/*
 * Marathon Spatial physics audio, C API.
 *
 * Modal bodies, contacts, debris and recorded layers, as mono streams for spatial sources. Register models, surfaces and
 * banks first, on one thread, before the first voice is pulled. After that, one control thread posts events and the
 * render thread pulls voices; neither locks or allocates, except ms_physics_fracture, which plans on the control thread.
 *
 *   ms_source_set_callback(engine, source, ms_physics_source_callback, ms_physics_voice(physics, slot));
 *
 * Nothing here opens an audio device.
 */
#ifndef MARATHON_SPATIAL_PHYSICS_H
#define MARATHON_SPATIAL_PHYSICS_H

#include "marathon_spatial.h"
#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct ms_physics ms_physics;
typedef struct ms_sample ms_sample;
typedef struct ms_layers ms_layers;
typedef uint32_t ms_modal_model;
typedef uint32_t ms_surface;
typedef uint32_t ms_debris_bank;

typedef struct ms_physics_desc {
    float sample_rate;
    uint32_t voices;            /* Modal voices, up to 64. */
    uint32_t contacts;          /* Sustained contacts at once, up to voices. */
    uint32_t modes_per_voice;
    uint32_t max_block_frames;  /* Largest pull, up to 8192. */
} ms_physics_desc;
void ms_physics_desc_default(ms_physics_desc* desc);
ms_status ms_physics_create(const ms_physics_desc* desc, ms_physics** physics);
void ms_physics_destroy(ms_physics* physics);

/* Plain data. Fill one yourself, or start from a preset and change any field. */
typedef struct ms_physics_material {
    float young_modulus_pa, density_kg_m3, poisson_ratio, loss_factor;
    float restitution, static_friction, dynamic_friction;
    float corner_wavelength_m, slide_fluctuation, roll_fluctuation, stick_stiffness_n_m;
    float panel_thickness_m;
} ms_physics_material;
uint32_t ms_physics_material_preset_count(void);
const char* ms_physics_material_preset_name(uint32_t index);
ms_status ms_physics_material_preset(const char* name, ms_physics_material* material);

/* Closed-form bar, plate or block, chosen from the box's proportions. wall_thickness_m above 0 makes a thin panel. */
ms_status ms_physics_model_box(ms_physics* physics, const ms_physics_material* material, ms_vec3 size_m, float wall_thickness_m,
                               ms_modal_model* model);
/* Finite-element analysis of a closed mesh, slow (seconds). Vertices are xyz triples in metres. */
ms_status ms_physics_model_mesh(ms_physics* physics, const ms_physics_material* material, const float* vertices, uint32_t vertex_count,
                                const uint32_t* indices, uint32_t index_count, float lattice_spacing_m, uint32_t max_modes,
                                ms_modal_model* model);
ms_status ms_physics_model_load(ms_physics* physics, const void* encoded, size_t size, ms_modal_model* model);
ms_status ms_physics_surface_create(ms_physics* physics, float corner_wavelength_m, uint32_t seed, ms_surface* surface);
ms_status ms_physics_debris_bank_create(ms_physics* physics, const ms_physics_material* material, float min_mass_kg, float max_mass_kg,
                                        uint32_t count, ms_debris_bank* bank);

/* J = (1 + e) m* v and the Hertz contact time. mass_b_kg may be INFINITY for a static body. */
ms_status ms_physics_solve_impact(const ms_physics_material* a, const ms_physics_material* b, float normal_speed_m_s, float mass_a_kg,
                                  float mass_b_kg, float radius_m, float restitution, float* impulse_ns, float* contact_seconds);

typedef struct ms_physics_impact_desc {
    uint64_t instance;          /* Nonzero key of the sounding object. */
    ms_modal_model model;
    ms_vec3 local_position;
    float impulse_ns, contact_seconds;
    float weight;               /* Above 0: distance gain at the listener, for voice competition. */
    float delay_seconds;
} ms_physics_impact_desc;
ms_status ms_physics_impact(ms_physics* physics, const ms_physics_impact_desc* impact);

typedef enum ms_contact_kind {
    MS_CONTACT_SLIDE = 0,
    MS_CONTACT_ROLL = 1,
    MS_CONTACT_SOURCE = 2       /* source and source_user give the force, scaled by force_scale_n. */
} ms_contact_kind;
typedef struct ms_physics_contact_desc {
    uint64_t contact;           /* Nonzero key of this touch. Two touches on a body are two contacts. */
    ms_modal_model model;
    ms_surface surface;
    ms_vec3 local_position;
    ms_contact_kind kind;
    float speed_m_s, normal_force_n;
    float dynamic_friction, static_friction, fluctuation, stick_stiffness_n_m;
    float radius_m, effective_modulus_pa, out_of_round;
    uint32_t flats;
    ms_layers* source;
    float force_scale_n;
    float ramp_seconds, weight;
} ms_physics_contact_desc;
void ms_physics_contact_desc_default(ms_physics_contact_desc* contact);
ms_status ms_physics_contact_update(ms_physics* physics, const ms_physics_contact_desc* contact);
ms_status ms_physics_contact_end(ms_physics* physics, uint64_t contact, float ramp_seconds);
ms_status ms_physics_stop(ms_physics* physics, uint64_t instance);

typedef struct ms_physics_fracture_desc {
    uint64_t id;                /* Nonzero. */
    ms_debris_bank bank;
    float parent_mass_kg, energy_j, height_m;
    uint32_t fragments, seed;
    const ms_physics_material* floor;
    float delay_seconds, weight;
} ms_physics_fracture_desc;
/* Plans the break and queues its hits. hits may be NULL. */
ms_status ms_physics_fracture(ms_physics* physics, const ms_physics_fracture_desc* fracture, uint32_t* hits);

/* The modal voice in a slot as a mono stream: pass this as user with ms_physics_source_callback. Pull sizes may differ
   between voices and blocks. Set doppler to 0 on these sources: the synth runs at one clock. */
void* ms_physics_voice(ms_physics* physics, uint32_t slot);
uint32_t ms_physics_source_callback(void* user, float* mono, uint32_t frames);
/* Runs the physics once per block: ms_engine_set_block_callback(engine, ms_physics_block_callback, physics). Without it
   nothing advances while no source is pulling, so a hit waits for a source to be playing. */
void ms_physics_block_callback(void* user, uint32_t frames, uint64_t frame);
/* The instance a slot sounds, so the host can place its source. 0 when idle. */
uint64_t ms_physics_voice_instance(const ms_physics* physics, uint32_t slot);
/* All voices as one mono mix, for hosts that do not place voices separately. */
void ms_physics_render_mix(ms_physics* physics, float* mono, uint32_t frames);

/* Bring your own audio: WAV, AIFF, CAF or anything else macOS decodes. */
ms_status ms_sample_load(const char* path, float sample_rate, ms_sample** sample);
void ms_sample_release(ms_sample* sample);
uint64_t ms_sample_frames(const ms_sample* sample);
int ms_sample_has_loop(const ms_sample* sample);

typedef enum ms_layer_kind {MS_LAYER_ONE_SHOT = 0, MS_LAYER_LOOP = 1, MS_LAYER_GRANULAR = 2} ms_layer_kind;
typedef enum ms_layer_target {
    MS_TARGET_GAIN = 0, MS_TARGET_PITCH = 1, MS_TARGET_LOW_PASS_HZ = 2, MS_TARGET_GRAIN_DENSITY = 3,
    MS_TARGET_GRAIN_SECONDS = 4, MS_TARGET_GRAIN_POSITION = 5, MS_TARGET_GRAIN_POSITION_JITTER = 6, MS_TARGET_GRAIN_PITCH_JITTER = 7
} ms_layer_target;
#define MS_PARAMETER_SPEED 0u
#define MS_PARAMETER_FORCE 1u
#define MS_PARAMETER_FIRST_USER 16u
#define MS_MAX_LAYER_MAPPINGS 8u
typedef struct ms_layer_mapping {
    uint32_t parameter;
    ms_layer_target target;
    float input_min, input_max, output_min, output_max, exponent;
    int logarithmic;
    float smoothing_seconds;
} ms_layer_mapping;
typedef struct ms_layer_desc {
    ms_layer_kind kind;
    const ms_sample* sample;
    float gain, pitch, low_pass_hz, crossfade_seconds;
    float grain_density, grain_seconds, grain_position, grain_position_jitter, grain_pitch_jitter;
    uint32_t mapping_count;
    ms_layer_mapping mappings[MS_MAX_LAYER_MAPPINGS];
} ms_layer_desc;
void ms_layer_desc_default(ms_layer_desc* layer);
ms_status ms_layers_create(const ms_layer_desc* layers, uint32_t count, float sample_rate, uint32_t seed, ms_layers** instrument);
void ms_layers_destroy(ms_layers* instrument);
void ms_layers_set_parameter(ms_layers* instrument, uint32_t parameter, float value);
void ms_layers_trigger(ms_layers* instrument);
void ms_layers_release(ms_layers* instrument, float seconds);
/* Pass this as user with ms_physics_source_callback. */
void* ms_layers_stream(ms_layers* instrument);

#ifdef __cplusplus
}
#endif
#endif
