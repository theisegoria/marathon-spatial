/* SPDX-FileCopyrightText: 2026 The Marathon Spatial Authors
 * SPDX-License-Identifier: Apache-2.0
 *
 * Writes interleaved float samples as 24-bit PCM, so the examples can show their result without an audio device.
 */
#ifndef EXAMPLES_WAV_WRITER_H
#define EXAMPLES_WAV_WRITER_H

#include <math.h>
#include <stdint.h>
#include <stdio.h>

static void wav_u32(FILE* file, uint32_t value) {
    unsigned char bytes[4] = {(unsigned char)value, (unsigned char)(value >> 8), (unsigned char)(value >> 16), (unsigned char)(value >> 24)};
    fwrite(bytes, 1, 4, file);
}

static void wav_u16(FILE* file, uint16_t value) {
    unsigned char bytes[2] = {(unsigned char)value, (unsigned char)(value >> 8)};
    fwrite(bytes, 1, 2, file);
}

/* Returns 0 on success. Samples outside [-1, 1] are clipped, so keep some headroom. */
static int wav_write(const char* path, const float* samples, uint32_t frames, uint16_t channels, uint32_t sample_rate) {
    FILE* file = fopen(path, "wb");
    if (!file) return -1;
    uint32_t data_bytes = frames * channels * 3u;
    fwrite("RIFF", 1, 4, file);
    wav_u32(file, 36u + data_bytes);
    fwrite("WAVEfmt ", 1, 8, file);
    wav_u32(file, 16);
    wav_u16(file, 1);
    wav_u16(file, channels);
    wav_u32(file, sample_rate);
    wav_u32(file, sample_rate * channels * 3u);
    wav_u16(file, (uint16_t)(channels * 3u));
    wav_u16(file, 24);
    fwrite("data", 1, 4, file);
    wav_u32(file, data_bytes);
    for (uint64_t index = 0; index < (uint64_t)frames * channels; ++index) {
        float clipped = samples[index] > 1.f ? 1.f : samples[index] < -1.f ? -1.f : samples[index];
        int32_t value = (int32_t)lroundf(clipped * 8388607.f);
        unsigned char bytes[3] = {(unsigned char)value, (unsigned char)(value >> 8), (unsigned char)(value >> 16)};
        fwrite(bytes, 1, 3, file);
    }
    return fclose(file) == 0 ? 0 : -1;
}

#endif
