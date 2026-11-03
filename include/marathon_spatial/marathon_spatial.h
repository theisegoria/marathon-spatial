/* SPDX-License-Identifier: Apache-2.0
 *
 * Marathon Spatial: spatial and physics-based audio for any headphones. macOS 27+, Apple silicon M3+.
 *
 * Threads
 *   Control side: every function here except the render entry points. Calls may come from any thread and are
 *   serialised by an internal lock that the render side never takes. They may allocate and may take time.
 *   Render side: ms_engine_render(), and any callback the engine invokes. It never locks, allocates or blocks.
 *   One thread at a time renders one engine.
 *
 * What is published when
 *   ms_*_set_* calls write a working copy. ms_engine_commit() publishes the whole working copy as one finished
 *   snapshot and hands the queued structural commands (play, stop, create, destroy) to the render side, which picks
 *   both up at the next block boundary. Nothing reaches the render side before a commit, except head orientation,
 *   which is a single atomic so head tracking keeps the lowest possible latency.
 *
 * Coordinates
 *   Right handed, metres, +X right, +Y up, -Z forward at the identity listener. Ambisonic content is ACN/SN3D (AmbiX),
 *   whose own X front, Y left, Z up frame the engine maps onto its axes.
 *
 * Determinism
 *   The same calls at the same rendered frame positions produce bit identical output, however the frames are chunked
 *   into ms_engine_render() calls. Live output and offline rendering use this one entry point.
 */
#ifndef MARATHON_SPATIAL_H
#define MARATHON_SPATIAL_H

#include <stddef.h>
#include <stdint.h>

#if defined(__cplusplus)
extern "C" {
#endif

#define MS_API_VERSION 1u
#define MS_MAX_AMBISONIC_ORDER 7u
#define MS_MAX_SENDS 2u
#define MS_MAX_REFLECTION_TAPS 24u
#define MS_MAX_SPEAKERS 32u

typedef enum ms_status {
    MS_OK = 0,
    MS_ERR_INVALID_ARGUMENT = -1,
    MS_ERR_OUT_OF_MEMORY = -2,
    MS_ERR_NOT_FOUND = -3,
    MS_ERR_BAD_FILE = -4,
    MS_ERR_LIMIT_REACHED = -5,
    MS_ERR_DEVICE = -6,
    MS_ERR_UNSUPPORTED = -7,
    MS_ERR_RATE_MISMATCH = -8
} ms_status;

typedef struct ms_vec3 { float x, y, z; } ms_vec3;
typedef struct ms_quat { float x, y, z, w; } ms_quat;

typedef struct ms_engine ms_engine;
typedef struct ms_hrtf ms_hrtf;

/* Handles. 0 is never valid. */
typedef uint32_t ms_source;
typedef uint32_t ms_buffer;
typedef uint32_t ms_bus;
typedef uint32_t ms_bed;
#define MS_BUS_MASTER ((ms_bus)1u)

/* Loudspeakers
 *
 * Without a layout the engine renders binaural stereo for headphones. With one it renders a channel per speaker: spatial sources pan
 * between the nearest speakers (vector base amplitude panning), and the ambisonic mix, which carries beds, reflections and the voices past
 * max_hrtf_voices, is decoded to the layout by AllRAD with max-rE weighting. The stereo master (pan sources, bus effects) plays as a phantom
 * pair at +-30 degrees. Head orientation is ignored: speakers stand still in the room.
 * Azimuth is in degrees anticlockwise from the front, so +30 is front left; elevation is up from ear height. LFE channels are left silent.
 * Presets follow ITU-R BS.2051 angles in the usual channel order (L R C LFE, then surrounds, then heights). */
typedef struct ms_speaker {
    float azimuth_degrees;
    float elevation_degrees;
    uint32_t lfe;  /* Nonzero for a low-frequency channel, which takes no panning. */
} ms_speaker;

typedef struct ms_speaker_layout {
    uint32_t count;
    ms_speaker speakers[MS_MAX_SPEAKERS];
} ms_speaker_layout;

typedef enum ms_speaker_preset {
    MS_SPEAKERS_MONO = 0,
    MS_SPEAKERS_STEREO = 1,   /* L R at +-30 */
    MS_SPEAKERS_QUAD = 2,     /* FL FR RL RR at +-45, +-135 */
    MS_SPEAKERS_5_1 = 3,      /* L R C LFE Ls Rs; surrounds at +-110 */
    MS_SPEAKERS_7_1 = 4,      /* L R C LFE Lss Rss Lrs Rrs; sides +-90, rears +-135 */
    MS_SPEAKERS_5_1_4 = 5,    /* 5.1, then heights at +-30 and +-110, 30 degrees up */
    MS_SPEAKERS_7_1_4 = 6     /* 7.1, then heights at +-45 and +-135, 30 degrees up */
} ms_speaker_preset;

ms_status ms_speaker_layout_preset(ms_speaker_preset preset, ms_speaker_layout* layout);
/* The panning gains the engine would give a direction, one per speaker, energy normalised. direction is in listener axes (+X right, +Y up,
   -Z forward). MS_ERR_INVALID_ARGUMENT for a layout the engine would refuse: none of it pannable, two speakers within a degree, too many. */
ms_status ms_speaker_layout_pan(const ms_speaker_layout* layout, ms_vec3 direction, float* gains);

/* Crosstalk cancellation
 *
 * Binaural sound on two loudspeakers. Each ear also hears the far speaker, and that crosstalk collapses the image to the line between the
 * speakers. The canceller feeds each speaker a filtered, inverted copy of the other channel that meets the crosstalk at the far ear, so the
 * binaural render reaches the ears nearly as on headphones and sounds can sit well outside the speakers. It works for a listener centred
 * between the speakers, facing them; off centre it fades towards plain stereo rather than failing badly, which the band limits and the
 * boost cap are for. Leave it off for headphones and for speaker layouts. */
typedef struct ms_crosstalk_params {
    uint32_t enabled;
    float speaker_span_degrees; /* Each speaker's angle off centre as the listener sees it: a laptop about 10, monitor speakers about 15,
                                   desktop speakers about 30. 2 to 60. */
    float head_radius_m;        /* The spherical head of the plant. 0.05 to 0.15. Default 0.0875. */
    float max_boost_db;         /* The most the canceller raises any frequency, 3 to 30. Lower is more robust and less coloured. Default 15. */
    float low_hz, high_hz;      /* The band where crosstalk is cancelled; plain stereo outside it. Default 200 and 7000. */
} ms_crosstalk_params;

void ms_crosstalk_params_default(ms_crosstalk_params* params);

/* Engine */

typedef struct ms_engine_desc {
    uint32_t api_version;       /* MS_API_VERSION */
    float sample_rate;          /* 8000 to 192000 */
    uint32_t block_frames;      /* 64, 128, 256 or 512 */
    uint32_t max_sources;       /* Sources that may exist. */
    uint32_t max_hrtf_voices;   /* Sources convolved directly. The rest go through the ambisonic mix. */
    uint32_t max_buses;
    uint32_t max_beds;
    uint32_t max_buffers;
    uint32_t ambisonic_order;   /* Shared mix order, 0 to MS_MAX_AMBISONIC_ORDER. 0 means no beds, no voice overflow. */
    float speed_of_sound;       /* m/s */
    float smoothing_seconds;    /* Pose smoothing time constant. 0 snaps at every block. */
    float air_temperature_c;    /* ISO 9613-1 air absorption */
    float air_humidity_percent;
    uint32_t limiter;           /* Nonzero: a master peak limiter at 0 dBFS. Off keeps measurements linear. */
    ms_speaker_layout speakers; /* count 0 for headphones; see ms_speaker_layout_preset. */
} ms_engine_desc;

void ms_engine_desc_default(ms_engine_desc* desc);
ms_status ms_engine_create(const ms_engine_desc* desc, ms_engine** engine);
void ms_engine_destroy(ms_engine* engine);
const char* ms_status_string(ms_status status);
const char* ms_version_string(void);

/* Publishes everything set since the last commit, and frees whatever the render side has finished with. All commands of
   one commit take effect in the same block, so sources played together with at_frame 0 start together. If the command
   queue lacks room for all of them they wait for a later commit; only a commit of more commands than the queue holds
   (8192) can be split. */
ms_status ms_engine_commit(ms_engine* engine);

/* Runs on the render side at the start of every block, whether or not anything is playing, before any source is pulled.
   Same rules as ms_source_callback. frame is the engine frame the block starts at. NULL removes it. The user data must
   stay alive as for ms_source_callback. */
typedef void (*ms_block_callback)(void* user, uint32_t frames, uint64_t frame);
ms_status ms_engine_set_block_callback(ms_engine* engine, ms_block_callback callback, void* user);

/* Renders exactly frames samples per channel into caller memory. Real-time safe. Returns the frame position after. */
uint64_t ms_engine_render(ms_engine* engine, float* left, float* right, uint32_t frames);
uint64_t ms_engine_render_interleaved(ms_engine* engine, float* stereo, uint32_t frames);
/* One channel per speaker of the engine's layout, in layout order. ms_engine_render on a speaker engine gives a stereo fold-down of the
   same render instead. Use one of the two per engine: a partial block left over from one is dropped by the other. Real-time safe. */
uint64_t ms_engine_render_speakers(ms_engine* engine, float* const* channels, uint32_t frames);
/* Crosstalk cancellation on the stereo master, for binaural sound on two speakers. Takes effect at the next commit, with a one-block
   crossfade. It adds latency_frames of delay while on. MS_ERR_INVALID_ARGUMENT for parameters out of range, MS_ERR_UNSUPPORTED on a
   speaker-layout engine. NULL turns it off. */
ms_status ms_engine_set_crosstalk(ms_engine* engine, const ms_crosstalk_params* params);
/* The delay the canceller adds while on, in frames; 0 when it is off. */
uint32_t ms_engine_crosstalk_latency_frames(const ms_engine* engine);
/* 2 for headphones, else the layout's speaker count. */
uint32_t ms_engine_output_channels(const ms_engine* engine);

/* The total latency of the direct path in frames, excluding any device latency. A set resampled to the engine's rate adds
   a further 8 frames' worth at 48 kHz (167 microseconds) to every direction alike. */
uint32_t ms_engine_direct_latency_frames(const ms_engine* engine);
float ms_engine_sample_rate(const ms_engine* engine);
uint32_t ms_engine_block_frames(const ms_engine* engine);

/* HRTF sets. A set is immutable, reference counted and may serve several engines. */

typedef struct ms_hrtf_info {
    char name[64];
    char licence[64];
    float sample_rate;
    float head_radius_m;
    float measured_distance_m;
    uint32_t directions;
    uint32_t taps;
} ms_hrtf_info;

typedef struct ms_hrtf_import_options {
    float sample_rate;          /* Target rate; 0 keeps the file's own. */
    uint32_t taps;              /* Minimum-phase filter length; 0 chooses 256 at 48 kHz, scaled with the rate. */
    uint32_t diffuse_field_eq;  /* Nonzero flattens the diffuse-field average, so sets sound alike in level and tone. */
    float head_radius_m;        /* 0 estimates from the set's interaural delays. */
} ms_hrtf_import_options;

void ms_hrtf_import_options_default(ms_hrtf_import_options* options);
/* Loads a cooked .mshrtf pack, resampling if sample_rate is nonzero and differs. */
ms_status ms_hrtf_load(const char* path, float sample_rate, ms_hrtf** hrtf);
/* Imports a user's AES69 SOFA file of type SimpleFreeFieldHRIR. Never call this from the render side. */
ms_status ms_hrtf_import_sofa(const char* path, const ms_hrtf_import_options* options, ms_hrtf** hrtf);
/* An analytic rigid sphere with ears at +-90 degrees. Needs no data; useful as a neutral fallback and a reference. */
ms_status ms_hrtf_create_sphere(float head_radius_m, float sample_rate, ms_hrtf** hrtf);
ms_status ms_hrtf_save(const ms_hrtf* hrtf, const char* path);
void ms_hrtf_get_info(const ms_hrtf* hrtf, ms_hrtf_info* info);
void ms_hrtf_retain(ms_hrtf* hrtf);
void ms_hrtf_release(ms_hrtf* hrtf);

/* Selects the set. Builds the ambisonic decoder, which can take a second, and crossfades over one block once committed. */
ms_status ms_engine_set_hrtf(ms_engine* engine, ms_hrtf* hrtf);

/* Listener */

ms_status ms_listener_set_pose(ms_engine* engine, ms_vec3 position, ms_vec3 velocity, ms_vec3 forward, ms_vec3 up);
/* Head orientation relative to the listener pose, from any tracker. Takes effect at the next block without a commit. */
void ms_head_set_orientation(ms_engine* engine, ms_quat orientation);

/* Sample buffers. Immutable once created; the engine owns a copy, resampled to its rate. */

ms_status ms_buffer_create(ms_engine* engine, const float* interleaved, uint64_t frames, uint32_t channels, float sample_rate,
                           ms_buffer* buffer);
/* Decodes any file Core Audio reads, on the calling thread. */
ms_status ms_buffer_load_file(ms_engine* engine, const char* path, ms_buffer* buffer);
ms_status ms_buffer_destroy(ms_engine* engine, ms_buffer buffer);
ms_status ms_buffer_get_info(ms_engine* engine, ms_buffer buffer, uint64_t* frames, uint32_t* channels);

/* Sources */

typedef enum ms_spatial_mode {
    MS_SPATIAL_HRTF = 0,        /* Direct convolution, falling back to the ambisonic mix past max_hrtf_voices. */
    MS_SPATIAL_AMBISONIC = 1,   /* Always the shared mix: cheap, for dense crowds and debris. */
    MS_SPATIAL_PAN = 2          /* Head locked equal-power stereo, using params.pan. For interface sounds. */
} ms_spatial_mode;

typedef enum ms_distance_model {
    MS_DISTANCE_NONE = 0,
    MS_DISTANCE_INVERSE = 1,    /* ref / (ref + rolloff (d - ref)) */
    MS_DISTANCE_LINEAR = 2,     /* 1 - rolloff (d - ref) / (max - ref) */
    MS_DISTANCE_EXPONENTIAL = 3 /* (d / ref) ^ -rolloff */
} ms_distance_model;

typedef struct ms_reflection_tap {
    float delay_seconds;        /* (reflected path - direct path) / speed of sound: how long after the direct sound it
                                   arrives, 0 to 0.17. Never the absolute path time. The tap reads the source's own signal,
                                   so this holds whether propagation_delay is on or off. */
    ms_vec3 direction;          /* World direction the reflection arrives from, seen from the listener. */
    float gain[3];              /* Linear, low (below 400 Hz), mid, high (above 2500 Hz). */
} ms_reflection_tap;

typedef struct ms_source_params {
    ms_vec3 position;
    ms_vec3 velocity;           /* m/s, for Doppler only. */
    ms_vec3 forward;            /* Directivity axis. */
    float gain;
    float pitch;                /* Playback ratio, on top of Doppler. 1/32 to 4. Content at another rate may use it. */
    uint32_t mode;              /* ms_spatial_mode */
    float pan;                  /* MS_SPATIAL_PAN: -1 left to +1 right. */
    uint32_t distance_model;    /* ms_distance_model */
    float ref_distance;
    float max_distance;
    float rolloff;
    float doppler;              /* 0 off, 1 physical. */
    float air_absorption;       /* 0 off, 1 ISO 9613-1. */
    float near_field;           /* 0 off, 1 rigid-sphere near-field model under the set's measured distance. */
    uint32_t propagation_delay; /* Nonzero delays the start by distance / speed of sound. */
    float directivity_alpha;    /* 0 omni, 0.5 cardioid, 1 dipole: (1 - a) + a cos(angle) */
    float directivity_sharpness;/* Exponent on the pattern. */
    float directivity_hf_db;    /* Extra high shelf loss at the rear, above 2 kHz. */
    float highpass_hz;          /* 0 off */
    float lowpass_hz;           /* 0 off */
    float eq_low_db;            /* Per-source transmission EQ for occlusion or material: below 400 Hz. */
    float eq_mid_db;            /* Around 1 kHz. */
    float eq_high_db;           /* Above 2.5 kHz. */
    ms_bus bus;                 /* 0 means master. */
    /* Sends hear the source with gain, distance and directivity, but before highpass, lowpass, EQ and air: those filter
       the path to the listener, and a source behind a door still excites its own room. */
    ms_bus send_bus[MS_MAX_SENDS];
    float send_gain[MS_MAX_SENDS];
} ms_source_params;

/* Pulls mono frames on the render side. Returns the frames written; fewer than asked ends the source. Must obey the
   render side rules. Destroying a source does not wait for the render side: keep user data alive until
   ms_engine_get_stats shows two more blocks rendered after the commit, or until the engine is destroyed. */
typedef uint32_t (*ms_source_callback)(void* user, float* mono, uint32_t frames);

void ms_source_params_default(ms_source_params* params);
/* Every float must be finite, including those of each ms_reflection_tap; anything else is MS_ERR_INVALID_ARGUMENT. */
ms_status ms_source_create(ms_engine* engine, ms_source* source);
ms_status ms_source_destroy(ms_engine* engine, ms_source source);
ms_status ms_source_set_params(ms_engine* engine, ms_source source, const ms_source_params* params);
ms_status ms_source_get_params(ms_engine* engine, ms_source source, ms_source_params* params);
ms_status ms_source_set_pose(ms_engine* engine, ms_source source, ms_vec3 position, ms_vec3 velocity);
ms_status ms_source_set_buffer(ms_engine* engine, ms_source source, ms_buffer buffer, uint32_t channel, int loop);
ms_status ms_source_set_callback(ms_engine* engine, ms_source source, ms_source_callback callback, void* user);
/* at_frame is an absolute engine frame; earlier frames, or 0, mean the next block. Sample accurate. */
ms_status ms_source_play(ms_engine* engine, ms_source source, uint64_t at_frame);
ms_status ms_source_stop(ms_engine* engine, ms_source source);
ms_status ms_source_seek(ms_engine* engine, ms_source source, double frame);
/* Early reflections for this source, replacing the previous list. Rendered through the ambisonic mix with its own
   band-split delay lines, so a ray tracer can publish many per source. A slot index keeps its path: delay and gain glide
   from the old value, so a path that persists must stay in its slot. Direct voices only. Published at the next commit. */
ms_status ms_source_set_reflections(ms_engine* engine, ms_source source, const ms_reflection_tap* taps, uint32_t count);

typedef struct ms_source_state {
    uint32_t playing;           /* Includes the filter tail after the content ends. */
    uint32_t direct;            /* Nonzero while convolved directly. */
    double position_frames;     /* Read position in the content. */
    float level;                /* Linear peak of the last block after distance and filters. */
} ms_source_state;

ms_status ms_source_get_state(ms_engine* engine, ms_source source, ms_source_state* state);

/* Ambisonic beds */

/* Pulls interleaved ambisonic frames on the render side; same rules as ms_source_callback. */
typedef uint32_t (*ms_bed_callback)(void* user, float* interleaved, uint32_t channels, uint32_t frames);

ms_status ms_bed_create(ms_engine* engine, ms_bed* bed);
ms_status ms_bed_destroy(ms_engine* engine, ms_bed bed);
/* The buffer holds (order + 1)^2 ACN/SN3D channels. Orders above the engine's are truncated. */
ms_status ms_bed_set_buffer(ms_engine* engine, ms_bed bed, ms_buffer buffer, int loop);
ms_status ms_bed_set_callback(ms_engine* engine, ms_bed bed, ms_bed_callback callback, void* user, uint32_t order);
ms_status ms_bed_set_gain(ms_engine* engine, ms_bed bed, float gain);
/* Scales the bed by this bus's gain and its parents', as for a source's ambisonic path. 0 means master. */
ms_status ms_bed_set_bus(ms_engine* engine, ms_bed bed, ms_bus bus);
/* The bed's orientation in the world. head_locked keeps it fixed to the head instead of the world. */
ms_status ms_bed_set_orientation(ms_engine* engine, ms_bed bed, ms_quat orientation, int head_locked);
ms_status ms_bed_play(ms_engine* engine, ms_bed bed, uint64_t at_frame);
ms_status ms_bed_stop(ms_engine* engine, ms_bed bed);

/* Buses. A bus mixes stereo into its parent; parents must already exist. MS_BUS_MASTER always does. Ambisonic content
   (overflow voices, reflections, beds) joins the master, scaled by the gains of its bus and that bus's parents, and passes
   through the master's own gain once, like everything else. */

/* In-place stereo processing on the render side. A reverb lives here. */
typedef void (*ms_bus_processor)(void* user, float* left, float* right, uint32_t frames, uint64_t frame);

ms_status ms_bus_create(ms_engine* engine, ms_bus parent, ms_bus* bus);
ms_status ms_bus_set_gain(ms_engine* engine, ms_bus bus, float gain);
ms_status ms_bus_set_filter(ms_engine* engine, ms_bus bus, float highpass_hz, float lowpass_hz);
ms_status ms_bus_set_processor(ms_engine* engine, ms_bus bus, ms_bus_processor processor, void* user);

/* Telemetry */

typedef struct ms_engine_stats {
    uint64_t frames_rendered;
    uint64_t blocks_rendered;
    uint32_t sources_playing;
    uint32_t voices_direct;
    uint32_t voices_ambisonic;
    uint32_t commands_pending;
    float peak[2];              /* Of the last block */
    float rms[2];
    uint64_t overs;             /* Output samples above full scale so far. */
    uint64_t faults;            /* Blocks of one source, bed or bus, or of the master as a last resort, replaced by silence
                                   because a sample was not finite. */
    double render_seconds_last; /* Wall clock. Telemetry only; it never affects the audio. */
    double render_seconds_max;
} ms_engine_stats;

/* Operations counted on the render side since creation, for a budget independent of the host. */
typedef struct ms_op_counts {
    uint64_t blocks;
    uint64_t fft_forward;       /* Real FFTs of 2 x block_frames */
    uint64_t fft_inverse;
    uint64_t spectrum_macs;     /* Complex multiply-accumulates of one bin */
    uint64_t real_macs;         /* Real multiply-accumulates: filters, interpolators, ramps, rotations */
    uint64_t filter_updates;    /* Filters resynthesised because a direction changed */
} ms_op_counts;

void ms_engine_get_stats(const ms_engine* engine, ms_engine_stats* stats);
void ms_engine_get_op_counts(const ms_engine* engine, ms_op_counts* counts);

#if defined(__cplusplus)
}
#endif
#endif
