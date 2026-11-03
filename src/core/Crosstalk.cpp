// SPDX-License-Identifier: Apache-2.0
#include "Crosstalk.hpp"

#include <algorithm>
#include <cmath>

#include "Fft.hpp"

namespace mspatial::core {
namespace {

constexpr double kPiD = 3.14159265358979323846;

// Brown and Duda's head shadow at the ear for a source theta radians off that ear's axis.
std::complex<double> shadow(double omega, double theta, double radius, double c) {
    constexpr double alphaMin = 0.1, thetaMin = 150.0 * kPiD / 180.0;
    const double alpha = (1.0 + alphaMin / 2.0) + (1.0 - alphaMin / 2.0) * std::cos(theta / thetaMin * kPiD);
    const double w0 = c / radius;
    return std::complex<double>(1.0, alpha * omega / (2.0 * w0)) / std::complex<double>(1.0, omega / (2.0 * w0));
}

// 1 inside [low, high], falling to 0 over an octave either side along a raised cosine.
double bandWeight(double hz, double low, double high) {
    const auto edge = [](double octaves) { return octaves >= 1.0 ? 0.0 : 0.5 * (1.0 + std::cos(kPiD * octaves)); };
    if (hz <= 0.0) return 0.0;
    if (hz < low) return edge(std::log2(low / hz));
    if (hz > high) return edge(std::log2(hz / high));
    return 1.0;
}

}  // namespace

std::complex<double> crosstalkPlant(const ms_crosstalk_params& p, double hz, double c) {
    const double span = p.speaker_span_degrees * kPiD / 180.0, a = p.head_radius_m, omega = 2.0 * kPiD * hz;
    const double delay = a / c * (span + std::sin(span));  // Woodworth: the far ear hears it this much later.
    const std::complex<double> ratio = shadow(omega, kPiD / 2.0 + span, a, c) / shadow(omega, kPiD / 2.0 - span, a, c);
    return ratio * std::polar(1.0, -omega * delay);
}

std::unique_ptr<CrosstalkFilter> designCrosstalk(const ms_crosstalk_params& p, float rate, std::uint32_t B, double c) {
    if (!(p.speaker_span_degrees >= 2.0f && p.speaker_span_degrees <= 60.0f) || !(p.head_radius_m >= 0.05f && p.head_radius_m <= 0.15f) ||
        !(p.max_boost_db >= 3.0f && p.max_boost_db <= 30.0f) || !(p.low_hz >= 20.0f) || !(p.high_hz > p.low_hz * 2.0f) || !(rate > 0.0f) || B == 0)
        return nullptr;
    auto f = std::make_unique<CrosstalkFilter>();
    f->taps = rate > 64000.0f ? 4096 : 2048;  // About 43 ms: fine enough in the bass that the band edge stays put.
    f->blockFrames = B;
    f->partitions = (f->taps + B - 1) / B;
    // The modelling delay. The regularised inverse rings on both sides of its main tap, longest in the bass, where the two ears' paths differ
    // least; an eighth of the filter, 5.3 ms at 48 kHz, holds that ringing (40 dB of cancellation in the band, within a decibel of plain
    // stereo outside it), where a sixteenth gives up about 7 dB.
    f->latencyFrames = f->taps / 8;
    // Design on a grid four times the filter. The band is applied to the plant, not to the canceller: it cancels the crosstalk w R, with w
    // the band weight, so outside the band it is plain stereo by construction and nothing has to be blended. Blending a finished inverse
    // with plain stereo instead crosses two responses of very different phase at the band edge, which rings far longer than the filter.
    // Each channel's inverse is Tikhonov's x* / (|x|^2 + beta), x = 1 + w R or 1 - w R, which keeps below the boost cap because
    // |x| / (|x|^2 + beta) peaks at 1 / (2 sqrt(beta)).
    const std::uint32_t N = 4 * f->taps, bins = N / 2;
    const double cap = std::pow(10.0, p.max_boost_db / 20.0), beta = 1.0 / (4.0 * cap * cap);
    std::vector<std::complex<double>> sum(bins + 1), difference(bins + 1);
    for (std::uint32_t k = 0; k <= bins; ++k) {
        const double hz = static_cast<double>(k) * rate / N;
        const std::complex<double> r = bandWeight(hz, p.low_hz, p.high_hz) * crosstalkPlant(p, hz, c);
        const std::complex<double> xs = 1.0 + r, xd = 1.0 - r;
        const std::complex<double> delay = std::polar(1.0, -2.0 * kPiD * k / N * f->latencyFrames);
        sum[k] = std::conj(xs) / (std::norm(xs) + beta) * delay;
        difference[k] = std::conj(xd) / (std::norm(xd) + beta) * delay;
    }
    // The impulse response from the Hermitian spectrum by the inverse FFT, which sums without the 1/N; then a half-Hann fade over the last
    // eighth of the taps kept.
    RealFft design(N);
    std::vector<float> re(bins), im(bins), time(N);
    const auto impulse = [&](const std::vector<std::complex<double>>& spectrum, std::vector<float>& out) {
        re[0] = static_cast<float>(spectrum[0].real());
        im[0] = static_cast<float>(spectrum[bins].real());
        for (std::uint32_t k = 1; k < bins; ++k) {
            re[k] = static_cast<float>(spectrum[k].real());
            im[k] = static_cast<float>(spectrum[k].imag());
        }
        design.inverse(re.data(), im.data(), time.data());
        out.assign(f->taps, 0.0f);
        const std::uint32_t fade = f->taps / 8;
        for (std::uint32_t n = 0; n < f->taps; ++n) {
            double g = 1.0;
            if (n >= f->taps - fade) g = 0.5 * (1.0 + std::cos(kPiD * (n - (f->taps - fade) + 1) / fade));
            out[n] = static_cast<float>(time[n] / N * g);
        }
    };
    impulse(sum, f->sum);
    impulse(difference, f->difference);
    // Partition spectra, as the engine's convolutions take them: each block of taps zero-padded to two blocks.
    RealFft fft(2 * B);
    f->re.reset(2ull * f->partitions * B);
    f->im.reset(2ull * f->partitions * B);
    std::vector<float> segment(2 * B);
    for (std::uint32_t channel = 0; channel < 2; ++channel) {
        const std::vector<float>& taps = channel == 0 ? f->sum : f->difference;
        for (std::uint32_t part = 0; part < f->partitions; ++part) {
            std::fill(segment.begin(), segment.end(), 0.0f);
            for (std::uint32_t t = 0; t < B && part * B + t < f->taps; ++t) segment[t] = taps[part * B + t];
            const std::size_t at = (static_cast<std::size_t>(channel) * f->partitions + part) * B;
            fft.forward(segment.data(), f->re.data() + at, f->im.data() + at);
        }
    }
    return f;
}

}  // namespace mspatial::core
