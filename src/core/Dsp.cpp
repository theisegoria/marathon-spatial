// SPDX-License-Identifier: Apache-2.0
#include "Dsp.hpp"

#include <algorithm>
#include <vector>

namespace mspatial::core {
namespace {

double besselI0(double x) {
    double sum = 1.0, term = 1.0;
    for (int k = 1; k < 64; ++k) {
        term *= (x / (2.0 * k)) * (x / (2.0 * k));
        sum += term;
        if (term < 1e-15 * sum) break;
    }
    return sum;
}

double kaiser(double t, double beta) {  // t in [-1, 1]
    const double a = 1.0 - t * t;
    return a <= 0.0 ? 0.0 : besselI0(beta * std::sqrt(a)) / besselI0(beta);
}

double sinc(double x) { return std::fabs(x) < 1e-12 ? 1.0 : std::sin(M_PI * x) / (M_PI * x); }

SincTable* buildSinc() {
    auto* t = new SincTable();
    constexpr double beta = 7.0;
    for (int p = 0; p <= SincTable::kPhases; ++p) {
        const double fraction = static_cast<double>(p) / SincTable::kPhases;
        double sum = 0.0;
        double taps[SincTable::kTaps];
        for (int i = 0; i < SincTable::kTaps; ++i) {
            const double x = (i - 3) - fraction;
            taps[i] = sinc(x) * kaiser(x / 4.0, beta);
            sum += taps[i];
        }
        for (int i = 0; i < SincTable::kTaps; ++i) t->taps[p][i] = static_cast<float>(taps[i] / sum);
    }
    // Exact at phase zero, so an unmoved source passes through untouched.
    for (int i = 0; i < SincTable::kTaps; ++i) {
        t->taps[0][i] = i == 3 ? 1.0f : 0.0f;
        t->taps[SincTable::kPhases][i] = i == 4 ? 1.0f : 0.0f;
    }
    for (int p = 0; p < SincTable::kPhases; ++p)
        for (int i = 0; i < SincTable::kTaps; ++i) t->delta[p][i] = t->taps[p + 1][i] - t->taps[p][i];
    return t;
}

Biquad normalise(double b0, double b1, double b2, double a0, double a1, double a2) noexcept {
    Biquad q;
    q.b0 = static_cast<float>(b0 / a0);
    q.b1 = static_cast<float>(b1 / a0);
    q.b2 = static_cast<float>(b2 / a0);
    q.a1 = static_cast<float>(a1 / a0);
    q.a2 = static_cast<float>(a2 / a0);
    return q;
}

}  // namespace

const SincTable* const gSincTable = buildSinc();

Biquad designLowpass(float hz, float q, float rate) noexcept {
    const double w = 2.0 * M_PI * std::min<double>(hz, rate * 0.49) / rate, c = std::cos(w), alpha = std::sin(w) / (2.0 * q);
    return normalise((1 - c) / 2, 1 - c, (1 - c) / 2, 1 + alpha, -2 * c, 1 - alpha);
}
Biquad designHighpass(float hz, float q, float rate) noexcept {
    const double w = 2.0 * M_PI * std::min<double>(hz, rate * 0.49) / rate, c = std::cos(w), alpha = std::sin(w) / (2.0 * q);
    return normalise((1 + c) / 2, -(1 + c), (1 + c) / 2, 1 + alpha, -2 * c, 1 - alpha);
}
Biquad designLowShelf(float hz, float db, float rate) noexcept {
    const double A = std::pow(10.0, db / 40.0), w = 2.0 * M_PI * std::min<double>(hz, rate * 0.49) / rate;
    const double c = std::cos(w), alpha = std::sin(w) * M_SQRT1_2, s = 2 * std::sqrt(A) * alpha;
    return normalise(A * ((A + 1) - (A - 1) * c + s), 2 * A * ((A - 1) - (A + 1) * c), A * ((A + 1) - (A - 1) * c - s),
                     (A + 1) + (A - 1) * c + s, -2 * ((A - 1) + (A + 1) * c), (A + 1) + (A - 1) * c - s);
}
Biquad designHighShelf(float hz, float db, float rate) noexcept {
    const double A = std::pow(10.0, db / 40.0), w = 2.0 * M_PI * std::min<double>(hz, rate * 0.49) / rate;
    const double c = std::cos(w), alpha = std::sin(w) * M_SQRT1_2, s = 2 * std::sqrt(A) * alpha;
    return normalise(A * ((A + 1) + (A - 1) * c + s), -2 * A * ((A - 1) + (A + 1) * c), A * ((A + 1) + (A - 1) * c - s),
                     (A + 1) - (A - 1) * c + s, 2 * ((A - 1) - (A + 1) * c), (A + 1) - (A - 1) * c - s);
}
Biquad designPeak(float hz, float q, float db, float rate) noexcept {
    const double A = std::pow(10.0, db / 40.0), w = 2.0 * M_PI * std::min<double>(hz, rate * 0.49) / rate;
    const double c = std::cos(w), alpha = std::sin(w) / (2.0 * q);
    return normalise(1 + alpha * A, -2 * c, 1 - alpha * A, 1 + alpha / A, -2 * c, 1 - alpha / A);
}
Biquad designFirstOrderLowShelf(float hz, float lowDb, float rate) noexcept {
    // Bilinear transform of H(s) = (s + g w) / (s + w), prewarped.
    const double g = std::pow(10.0, lowDb / 20.0), k = std::tan(M_PI * std::min<double>(hz, rate * 0.49) / rate);
    const double a0 = 1 + k;
    Biquad q;
    q.b0 = static_cast<float>((1 + g * k) / a0);
    q.b1 = static_cast<float>((g * k - 1) / a0);
    q.b2 = 0;
    q.a1 = static_cast<float>((k - 1) / a0);
    q.a2 = 0;
    return q;
}

double airAttenuationDbPerMetre(double f, double celsius, double humidity, double kpa) noexcept {
    const double T = celsius + 273.15, T0 = 293.15, T01 = 273.16, pr = kpa / 101.325;
    const double psat = std::pow(10.0, -6.8346 * std::pow(T01 / T, 1.261) + 4.6151);
    const double h = humidity * psat / pr;
    const double frO = pr * (24.0 + 4.04e4 * h * (0.02 + h) / (0.391 + h));
    const double frN = pr * std::pow(T / T0, -0.5) * (9.0 + 280.0 * h * std::exp(-4.170 * (std::pow(T / T0, -1.0 / 3.0) - 1.0)));
    const double f2 = f * f;
    return 8.686 * f2 *
           (1.84e-11 / pr * std::sqrt(T / T0) +
            std::pow(T / T0, -2.5) * (0.01275 * std::exp(-2239.1 / T) / (frO + f2 / frO) +
                                      0.1068 * std::exp(-3352.0 / T) / (frN + f2 / frN)));
}

void resampleSinc(const float* input, std::uint64_t inputFrames, std::uint32_t stride, double inRate, double outRate,
                  float* output, std::uint64_t outputFrames, std::uint32_t leadFrames) noexcept {
    constexpr int halfTaps = 24;
    constexpr double beta = 10.0;
    const double cutoff = std::min(1.0, outRate / inRate) * 0.97;
    const double step = inRate / outRate;
    const double reach = halfTaps / cutoff;
    for (std::uint64_t n = 0; n < outputFrames; ++n) {
        const double centre = (static_cast<double>(n) - leadFrames) * step;
        const auto first = static_cast<std::int64_t>(std::ceil(centre - reach));
        const auto last = static_cast<std::int64_t>(std::floor(centre + reach));
        double sum = 0.0;
        for (std::int64_t i = std::max<std::int64_t>(first, 0); i <= last && i < static_cast<std::int64_t>(inputFrames); ++i) {
            const double x = static_cast<double>(i) - centre;
            sum += input[static_cast<std::uint64_t>(i) * stride] * sinc(x * cutoff) * cutoff * kaiser(x / reach, beta);
        }
        output[n] = static_cast<float>(sum);
    }
}

}  // namespace mspatial::core
