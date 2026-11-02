/* SPDX-FileCopyrightText: 2026 The Marathon Spatial Authors
 * SPDX-License-Identifier: Apache-2.0
 *
 * Places your own sound in a small room and renders a binaural WAV file offline. The input is any file Core Audio
 * reads. Pass a SOFA file of type SimpleFreeFieldHRIR to listen through your own HRTF; otherwise the analytic sphere
 * is used. The early reflections here are written by hand from the image-source method, to show what a ray tracer
 * publishes each frame. Nothing opens an audio device.
 *
 *     bring_your_own_audio INPUT OUTPUT.wav [HRTF.sofa]
 */
#include "marathon_spatial/marathon_spatial.h"
#include "wav_writer.h"

#include <stdio.h>
#include <stdlib.h>

#define SAMPLE_RATE 48000
#define BLOCK 256
#define TAIL_SECONDS 1

#define CHECK(call)                                                                  \
    do {                                                                             \
        ms_status status_ = (call);                                                  \
        if (status_ != MS_OK) {                                                      \
            fprintf(stderr, "%s: %s\n", #call, ms_status_string(status_));           \
            return 1;                                                                \
        }                                                                            \
    } while (0)

static float length(ms_vec3 v) { return sqrtf(v.x * v.x + v.y * v.y + v.z * v.z); }

/* First-order image of the source in a plane through `point` with unit normal `normal`. */
static ms_reflection_tap image_tap(ms_vec3 source, ms_vec3 point, ms_vec3 normal, float reflectance_low, float reflectance_high) {
    const float distance_to_plane = (source.x - point.x) * normal.x + (source.y - point.y) * normal.y + (source.z - point.z) * normal.z;
    const ms_vec3 image = {source.x - 2.f * distance_to_plane * normal.x, source.y - 2.f * distance_to_plane * normal.y,
                           source.z - 2.f * distance_to_plane * normal.z};
    const float path = length(image);
    const float mid = 0.5f * (reflectance_low + reflectance_high);
    ms_reflection_tap tap;
    tap.delay_seconds = (path - length(source)) / 343.f;
    tap.direction = (ms_vec3){image.x / path, image.y / path, image.z / path};
    tap.gain[0] = reflectance_low / path;
    tap.gain[1] = mid / path;
    tap.gain[2] = reflectance_high / path;
    return tap;
}

int main(int argument_count, char** arguments) {
    if (argument_count != 3 && argument_count != 4) {
        fprintf(stderr, "usage: %s INPUT OUTPUT.wav [HRTF.sofa]\n", arguments[0]);
        return 2;
    }

    ms_engine_desc desc;
    ms_engine_desc_default(&desc);
    desc.sample_rate = SAMPLE_RATE;
    desc.block_frames = BLOCK;
    ms_engine* engine = NULL;
    CHECK(ms_engine_create(&desc, &engine));

    ms_hrtf* hrtf = NULL;
    if (argument_count == 4) {
        ms_hrtf_import_options options;
        ms_hrtf_import_options_default(&options);
        options.sample_rate = SAMPLE_RATE;
        CHECK(ms_hrtf_import_sofa(arguments[3], &options, &hrtf));
    } else {
        CHECK(ms_hrtf_create_sphere(0.0875f, SAMPLE_RATE, &hrtf));
    }
    CHECK(ms_engine_set_hrtf(engine, hrtf));
    ms_hrtf_release(hrtf);

    ms_buffer buffer = 0;
    CHECK(ms_buffer_load_file(engine, arguments[1], &buffer));
    uint64_t content_frames = 0;
    uint32_t channels = 0;
    CHECK(ms_buffer_get_info(engine, buffer, &content_frames, &channels));

    /* Ahead and to the left. Multichannel input uses its first channel, because a point source is mono. */
    const ms_vec3 position = {-1.2f, 0.f, -1.8f};
    ms_source source = 0;
    CHECK(ms_source_create(engine, &source));
    CHECK(ms_source_set_buffer(engine, source, buffer, 0, 0));
    ms_source_params params;
    ms_source_params_default(&params);
    params.position = position;
    params.air_absorption = 1.f;
    CHECK(ms_source_set_params(engine, source, &params));

    /* A 6 x 3 x 5 m room around the listener: plaster walls, a carpeted floor. Slots stay stable, as paths would. */
    const ms_reflection_tap taps[] = {
        image_tap(position, (ms_vec3){3.f, 0.f, 0.f}, (ms_vec3){-1.f, 0.f, 0.f}, 0.9f, 0.8f),
        image_tap(position, (ms_vec3){-3.f, 0.f, 0.f}, (ms_vec3){1.f, 0.f, 0.f}, 0.9f, 0.8f),
        image_tap(position, (ms_vec3){0.f, 0.f, -2.5f}, (ms_vec3){0.f, 0.f, 1.f}, 0.9f, 0.8f),
        image_tap(position, (ms_vec3){0.f, 0.f, 2.5f}, (ms_vec3){0.f, 0.f, -1.f}, 0.9f, 0.8f),
        image_tap(position, (ms_vec3){0.f, -1.5f, 0.f}, (ms_vec3){0.f, 1.f, 0.f}, 0.8f, 0.3f),
        image_tap(position, (ms_vec3){0.f, 1.5f, 0.f}, (ms_vec3){0.f, -1.f, 0.f}, 0.9f, 0.8f),
    };
    CHECK(ms_source_set_reflections(engine, source, taps, (uint32_t)(sizeof taps / sizeof taps[0])));
    CHECK(ms_source_play(engine, source, 0));
    CHECK(ms_engine_commit(engine));

    if (content_frames > UINT32_MAX / 4 - SAMPLE_RATE * TAIL_SECONDS) {
        fprintf(stderr, "input too long for this example\n");
        return 1;
    }
    const uint32_t frames = (uint32_t)content_frames + SAMPLE_RATE * TAIL_SECONDS;
    float* output = malloc((size_t)frames * 2 * sizeof(float));
    if (!output) return 1;
    ms_engine_render_interleaved(engine, output, frames);
    ms_engine_destroy(engine);

    int result = wav_write(arguments[2], output, frames, 2, SAMPLE_RATE);
    free(output);
    if (result != 0) {
        fprintf(stderr, "could not write %s\n", arguments[2]);
        return 1;
    }
    printf("rendered %u frames, wrote %s\n", frames, arguments[2]);
    return 0;
}
