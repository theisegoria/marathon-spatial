// SPDX-License-Identifier: Apache-2.0
// Sample-level building blocks of the render graph.
#pragma once
#include <cmath>
#include <cstdint>

#include "Common.hpp"

namespace mspatial::core {

// Eight-tap Kaiser-windowed sinc over 1024 phases. Taps sit at offsets -3..+4 around floor(position).
struct SincTable {
    static constexpr int kTaps = 8, kPhases = 1024;
    alignas(64) float taps[kPhases + 1][kTaps];
    float delta[kPhases][kTaps];
};
extern const SincTable* const gSincTable;  // Built at load, so the hot path needs no guard.

inline float sincSample(const float* centre, float fraction) noexcept {  // centre[-3..+4] readable
    const SincTable& t = *gSincTable;
    const float scaled = fraction * SincTable::kPhases;
    const int phase = static_cast<int>(scaled);
    const float within = scaled - static_cast<float>(phase);
    const float* k = t.taps[phase];
    const float* d = t.delta[phase];
    float sum = 0.0f;
    for (int i = 0; i < SincTable::kTaps; ++i) sum += (k[i] + within * d[i]) * centre[i - 3];
    return sum;
}

// Direct form II transposed biquad. Stepping coefficients leaves the state valid, which is why this form is used.
struct Biquad {
    float b0{1}, b1{}, b2{}, a1{}, a2{}, z1{}, z2{};
    bool identityCoefficients() const noexcept { return b0 == 1 && b1 == 0 && b2 == 0 && a1 == 0 && a2 == 0; }
    bool silentState() const noexcept { return std::fabs(z1) + std::fabs(z2) < 1e-12f; }
    void process(float* samples, std::uint32_t n) noexcept {
        float s1 = z1, s2 = z2;
        for (std::uint32_t i = 0; i < n; ++i) {
            const float x = samples[i];
            const float y = b0 * x + s1;
            s1 = b1 * x - a1 * y + s2;
            s2 = b2 * x - a2 * y;
            samples[i] = y;
        }
        z1 = s1;
        z2 = s2;
    }
    void setCoefficients(const Biquad& c) noexcept { b0 = c.b0; b1 = c.b1; b2 = c.b2; a1 = c.a1; a2 = c.a2; }
    void clear() noexcept { z1 = z2 = 0.0f; }
};

// Audio EQ Cookbook designs. Frequencies are normalised by the sample rate at the call.
Biquad designLowpass(float hz, float q, float rate) noexcept;
Biquad designHighpass(float hz, float q, float rate) noexcept;
Biquad designLowShelf(float hz, float db, float rate) noexcept;
Biquad designHighShelf(float hz, float db, float rate) noexcept;
Biquad designPeak(float hz, float q, float db, float rate) noexcept;
// A first-order low shelf as a biquad: gain lowDb at DC, 0 dB at high frequency.
Biquad designFirstOrderLowShelf(float hz, float lowDb, float rate) noexcept;

// ISO 9613-1 attenuation of sound in air, dB per metre.
double airAttenuationDbPerMetre(double hz, double celsius, double humidityPercent, double pressureKpa = 101.325) noexcept;

// Band-limited change of rate with a Kaiser-windowed sinc. Off the render side only. leadFrames of output come before the
// input's first sample, which keeps the part of a sudden onset that a band-limited signal has before it.
void resampleSinc(const float* input, std::uint64_t inputFrames, std::uint32_t stride, double inRate, double outRate,
                  float* output, std::uint64_t outputFrames, std::uint32_t leadFrames = 0) noexcept;
inline std::uint64_t resampledLength(std::uint64_t frames, double inRate, double outRate) noexcept {
    return static_cast<std::uint64_t>(std::ceil(static_cast<double>(frames) * outRate / inRate));
}

}  // namespace mspatial::core
