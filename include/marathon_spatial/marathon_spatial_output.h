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

/* Binds an engine to a device. device_uid is a Core Audio device UID, or NULL for the default output device.
   Fails with MS_ERR_RATE_MISMATCH when the device's rate is not the engine's. Nothing runs until ms_output_start(). */
ms_status ms_output_create(ms_engine* engine, const char* device_uid, ms_output** output);
ms_status ms_output_start(ms_output* output);
ms_status ms_output_stop(ms_output* output);
void ms_output_destroy(ms_output* output);

typedef struct ms_output_stats {
    uint64_t callbacks;
    uint64_t overruns;      /* Cycles the HAL reported as overloaded. */
    uint32_t buffer_frames;
    uint32_t workgroup_joined;
    uint32_t device_alive;
    uint32_t latency_frames; /* Device, stream and safety offset latency. */
} ms_output_stats;

void ms_output_get_stats(const ms_output* output, ms_output_stats* stats);

#if defined(__cplusplus)
}
#endif
#endif
