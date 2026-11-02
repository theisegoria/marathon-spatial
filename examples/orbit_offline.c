/* SPDX-FileCopyrightText: 2026 The Marathon Spatial Authors
 * SPDX-License-Identifier: Apache-2.0
 *
 * The smallest complete use of the library: one source orbiting the listener's head, rendered offline to a binaural
 * WAV file. The HRTF is the built-in analytic sphere, so this needs no data files. Nothing opens an audio device.
 *
 *     orbit_offline OUTPUT.wav
 */
#include "marathon_spatial/marathon_spatial.h"
#include "wav_writer.h"

#include <stdio.h>
#include <stdlib.h>

#define SAMPLE_RATE 48000
#define BLOCK 256
#define SECONDS 8
#define TWO_PI 6.28318530717958647692f

#define CHECK(call)                                                                  \
    do {                                                                             \
        ms_status status_ = (call);                                                  \
        if (status_ != MS_OK) {                                                      \
            fprintf(stderr, "%s: %s\n", #call, ms_status_string(status_));           \
            return 1;                                                                \
        }                                                                            \
    } while (0)

int main(int argument_count, char** arguments) {
    if (argument_count != 2) {
        fprintf(stderr, "usage: %s OUTPUT.wav\n", arguments[0]);
        return 2;
    }

    ms_engine_desc desc;
    ms_engine_desc_default(&desc);
    desc.sample_rate = SAMPLE_RATE;
    desc.block_frames = BLOCK;
    ms_engine* engine = NULL;
    CHECK(ms_engine_create(&desc, &engine));

    ms_hrtf* hrtf = NULL;
    CHECK(ms_hrtf_create_sphere(0.0875f, SAMPLE_RATE, &hrtf));
    CHECK(ms_engine_set_hrtf(engine, hrtf));
    ms_hrtf_release(hrtf); /* The engine keeps its own reference. */

    /* Broadband noise gives the HRTF energy to act on at every frequency. */
    const uint32_t content_frames = SAMPLE_RATE;
    float* noise = malloc(content_frames * sizeof(float));
    if (!noise) return 1;
    for (uint32_t index = 0, state = 22222; index < content_frames; ++index) {
        state = state * 1664525u + 1013904223u;
        noise[index] = 0.25f * ((float)(state >> 8) / 8388608.f - 1.f);
    }
    ms_buffer buffer = 0;
    CHECK(ms_buffer_create(engine, noise, content_frames, 1, SAMPLE_RATE, &buffer));
    free(noise);

    ms_source source = 0;
    CHECK(ms_source_create(engine, &source));
    CHECK(ms_source_set_buffer(engine, source, buffer, 0, 1));
    ms_source_params params;
    ms_source_params_default(&params);
    params.mode = MS_SPATIAL_HRTF;
    params.near_field = 1.f;
    CHECK(ms_source_set_params(engine, source, &params));
    CHECK(ms_source_play(engine, source, 0));

    const uint32_t frames = SAMPLE_RATE * SECONDS;
    float* output = malloc((size_t)frames * 2 * sizeof(float));
    if (!output) return 1;
    const float radius = 1.5f;
    const float angular_speed = TWO_PI / 4.f; /* One orbit every four seconds. */
    for (uint32_t frame = 0; frame < frames; frame += BLOCK) {
        const float time = (float)frame / SAMPLE_RATE;
        const float angle = angular_speed * time;
        /* -Z is forward, so this starts in front and passes the right ear first. */
        ms_vec3 position = {radius * sinf(angle), 0.f, -radius * cosf(angle)};
        ms_vec3 velocity = {radius * angular_speed * cosf(angle), 0.f, radius * angular_speed * sinf(angle)};
        CHECK(ms_source_set_pose(engine, source, position, velocity));
        CHECK(ms_engine_commit(engine));
        const uint32_t count = frames - frame < BLOCK ? frames - frame : BLOCK;
        ms_engine_render_interleaved(engine, output + (size_t)frame * 2, count);
    }

    ms_engine_stats stats;
    ms_engine_get_stats(engine, &stats);
    ms_engine_destroy(engine);
    int result = wav_write(arguments[1], output, frames, 2, SAMPLE_RATE);
    free(output);
    if (result != 0) {
        fprintf(stderr, "could not write %s\n", arguments[1]);
        return 1;
    }
    printf("rendered %llu frames, %llu samples over full scale, wrote %s\n", (unsigned long long)stats.frames_rendered,
           (unsigned long long)stats.overs, arguments[1]);
    return 0;
}
