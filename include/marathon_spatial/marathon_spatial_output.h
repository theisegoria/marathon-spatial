/* SPDX-License-Identifier: Apache-2.0
 *
 * Live output. This is the only part of the library that touches an audio device, and it lives in its own library
 * target so offline users, and every test, link code that cannot open one.
 */
#ifndef MARATHON_SPATIAL_OUTPUT_H
#define MARATHON_SPATIAL_OUTPUT_H

#include "marathon_spatial.h"

#if defined(__cplusplus)
extern "C" {
#endif

typedef struct ms_output ms_output;

/* The nominal rate of the default output device, so an engine can be created to match it. */
ms_status ms_output_default_sample_rate(float* sample_rate);

/* Where the speakers of an output device stand, in its channel order, from the layout the device reports (its Audio MIDI Setup speaker
   configuration). Create an engine with this layout and the output plays one engine channel to each device channel. Channels with no
   position (unknown or discrete) come back as lfe and stay silent. device_uid as for ms_output_create. MS_ERR_NOT_FOUND when the device
   reports no layout it can place; a plain stereo device reports L and R. */
ms_status ms_output_device_layout(const char* device_uid, ms_speaker_layout* layout);

/* Binds an engine to a device. device_uid is a Core Audio device UID, or NULL for the default output device.
   Fails with MS_ERR_RATE_MISMATCH when the device's rate is not the engine's. Nothing runs until ms_output_start().
   A speaker engine plays engine channel c on device channel c when the device has a channel for each; otherwise, like a headphone engine,
   it plays stereo on the first two (a speaker engine folded down), and a mono device takes the mid signal. */
ms_status ms_output_create(ms_engine* engine, const char* device_uid, ms_output** output);
ms_status ms_output_start(ms_output* output);
ms_status ms_output_stop(ms_output* output);
void ms_output_destroy(ms_output* output);

/* Position tracking by the front camera, for steering crosstalk cancellation to where the listener sits. Apple's Vision framework finds
   the box round each person's upper body; the largest is followed, whoever else is in view, and the canceller keeps running. Only the
   position leaves the tracker: it never estimates which way anyone faces, and no image is kept. The first start asks the user for camera
   access, which the host app's Info.plist must explain under NSCameraUsageDescription. */
typedef struct ms_head_tracker ms_head_tracker;

typedef struct ms_head_tracker_config {
    float horizontal_fov_degrees;  /* The camera's horizontal field of view; macOS does not report it. 0 takes 70, about a Mac's built-in camera. */
    float body_width_m;            /* The real width of an upper body. 0 takes 0.42, about the shoulders. */
    float frames_per_second;       /* Frames looked at each second. 0 takes 15. */
    float smoothing;               /* How much of each new reading is taken, 0 to 1. 0 takes 0.35. */
} ms_head_tracker_config;

/* MS_ERR_DEVICE when the Mac has no camera. Nothing runs until ms_head_tracker_start. config may be NULL. */
ms_status ms_head_tracker_create(const ms_head_tracker_config* config, ms_head_tracker** tracker);
/* Starts the camera, asking for access the first time. MS_ERR_UNSUPPORTED when access has been refused. */
ms_status ms_head_tracker_start(ms_head_tracker* tracker);
void ms_head_tracker_stop(ms_head_tracker* tracker);
void ms_head_tracker_destroy(ms_head_tracker* tracker);
/* The latest smoothed position. valid is 0 until someone is seen, and again once no one has been for half a second. Any thread. */
void ms_head_tracker_get(const ms_head_tracker* tracker, ms_head_position* head);

typedef struct ms_output_stats {
    uint64_t callbacks;
    uint64_t overruns;      /* Cycles the HAL reported as overloaded. */
    uint32_t buffer_frames;
    uint32_t workgroup_joined;
    uint32_t device_alive;
    uint32_t latency_frames; /* Device, stream and safety offset latency. */
    uint32_t device_channels;  /* Output channels the device has. */
    uint32_t engine_channels;  /* Channels the engine renders to it: its speakers, or 2 when it plays stereo. */
} ms_output_stats;

void ms_output_get_stats(const ms_output* output, ms_output_stats* stats);

#if defined(__cplusplus)
}
#endif
#endif
