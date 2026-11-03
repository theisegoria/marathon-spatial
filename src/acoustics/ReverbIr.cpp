#include "acoustics/ReverbIr.hpp"

#include "core/Fft.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>

namespace mspatial::acoustics {
namespace {

constexpr float kPi = std::numbers::pi_v<float>;
constexpr float kUniformPrior = 0.1f;  // Share of each bin's energy spread evenly over directions, against noise in sparse bins.

// Power-complementary octave weights: squares sum to one at every frequency, so the bands of a white noise add back to white.
float bandWeight(std::uint32_t band, float frequency) {
    const float lowest = kBandCenterHz[0], highest = kBandCenterHz[kBandCount - 1];
    if (frequency <= lowest) return band == 0 ? 1.f : 0.f;
    if (frequency >= highest) return band == kBandCount - 1 ? 1.f : 0.f;
    const float octaves = std::log2(frequency / lowest);
    const auto below = std::min<std::uint32_t>(std::uint32_t(octaves), kBandCount - 2);
    const float x = octaves - float(below);
    if (band == below) return std::cos(0.5f * kPi * x);
    if (band == below + 1) return std::sin(0.5f * kPi * x);
    return 0.f;
}

float smoothstep(float edge0, float edge1, float x) {
    if (edge1 <= edge0) return x >= edge1 ? 1.f : 0.f;
    const float t = std::clamp((x - edge0) / (edge1 - edge0), 0.f, 1.f);
    return t * t * (3.f - 2.f * t);
}

}  // namespace

void ReverbIrSynth::prepare(float sampleRate, std::uint32_t maxSamples) {
    sampleRate_ = sampleRate;
    std::uint32_t size = 1024;
    while (size < maxSamples) size *= 2;
    // The bank is periodic in `size`: reading it at any circular offset never meets a seam.
    length_ = size;
    core::RealFft fft(size);
    std::vector<float> white(size), re(size / 2), im(size / 2);
    std::uint32_t state = 0x6d2b79f5u;
    auto uniform = [&state] {
        state ^= state << 13;
        state ^= state >> 17;
        state ^= state << 5;
        return (float(state >> 8) + 0.5f) * (1.f / 16777216.f);
    };
    noise_.assign(std::size_t(kBandCount) * size, 0.f);
    std::vector<float> whiteRe(size / 2), whiteIm(size / 2);
    for (std::uint32_t b = 0; b < kBandCount; ++b) {
        // Each band filters its own white noise. Bands cut from one noise would correlate where neighbours overlap, and their powers would
        // no longer add. Each noise has exactly unit mean square; the band keeps its share of that power, not a unit of its own: an impulse's
        // energy is shared across the bands, and the field's per-band values are each relative to the whole of the direct sound.
        for (std::uint32_t i = 0; i < size; i += 2) {
            const float r = std::sqrt(-2.f * std::log(uniform())), phi = 2.f * kPi * uniform();
            white[i] = r * std::cos(phi);
            white[i + 1] = r * std::sin(phi);
        }
        float whitePower = 0.f;
        vDSP_measqv(white.data(), 1, &whitePower, size);
        const float whiteScale = 1.f / std::sqrt(whitePower);
        vDSP_vsmul(white.data(), 1, &whiteScale, white.data(), 1, size);
        fft.forward(white.data(), whiteRe.data(), whiteIm.data());
        for (std::uint32_t k = 0; k < size / 2; ++k) {
            const float w = bandWeight(b, float(k) * sampleRate / float(size));
            re[k] = whiteRe[k] * w;
            im[k] = whiteIm[k] * w;
        }
        re[0] = whiteRe[0] * bandWeight(b, 0.f);
        im[0] = whiteIm[0] * bandWeight(b, 0.5f * sampleRate);
        float* out = noise_.data() + std::size_t(b) * size;
        fft.inverse(re.data(), im.data(), out);
        const float scale = 1.f / float(2 * size);  // forward then inverse gains 2N.
        vDSP_vsmul(out, 1, &scale, out, 1, size);
    }
    envelope_.assign(size, 0.f);
    signal_.assign(size, 0.f);
}

void ReverbIrSynth::smooth(const EnergyField& field) {
    smoothed_.bins = field.bins;
    smoothed_.binSeconds = field.binSeconds;
    smoothed_.room = field.room;
    smoothed_.paths = field.paths;
    smoothed_.energy.assign(field.energy.size(), 0.f);
    totals_.assign(field.bins, 0.f);
    // Double sums: a late bin is a millionth of an early one, which a float running sum would round away.
    std::vector<double> prefix(field.bins + 1);
    std::vector<double> directionPrefix(std::size_t(kFieldDirections) * (field.bins + 1));
    for (std::uint32_t b = 0; b < kBandCount; ++b) {
        prefix[0] = 0.0;
        for (std::uint32_t i = 0; i < field.bins; ++i) prefix[i + 1] = prefix[i] + field.band(b, i);
        for (std::uint32_t d = 0; d < kFieldDirections; ++d) {
            double* p = directionPrefix.data() + std::size_t(d) * (field.bins + 1);
            p[0] = 0.0;
            for (std::uint32_t i = 0; i < field.bins; ++i) p[i + 1] = p[i] + field.at(d, b, i);
        }
        for (std::uint32_t i = 0; i < field.bins; ++i) {
            // About a sixth of the elapsed time either side: tight where the early field has detail, wide where it is sparse and smooth.
            // The window stays centred, so it shrinks near the end.
            const std::uint32_t half = std::min({i / 6, 20u, field.bins - 1 - i});
            const std::uint32_t lo = i - half, hi = i + half + 1;
            const float count = float(hi - lo);
            const float mean = float((prefix[hi] - prefix[lo]) / count);
            float total = mean;
            // A plain mean over a decaying field overstates its centre (by sinh(x) / x for an exponential). The two halves give the local
            // decay, which corrects the mean back to the centre exactly for an exponential.
            if (half >= 2) {
                const float early = float((prefix[i] - prefix[lo]) / half), late = float((prefix[hi] - prefix[i + 1]) / half);
                if (early > 0.f && late > early * 1e-6f && late < early) {
                    const float step = std::log(early / late) / float(half + 1);  // Decay per bin, in nepers of energy.
                    float sum = 0.f;
                    for (int k = -int(half); k <= int(half); ++k) sum += std::exp(-step * float(k));
                    total *= count / sum;
                }
            }
            // Directions share the corrected total in the proportions the plain window holds.
            const float prior = kUniformPrior * mean / float(kFieldDirections);
            const float denominator = mean + float(kFieldDirections) * prior;
            for (std::uint32_t d = 0; d < kFieldDirections; ++d) {
                const double* p = directionPrefix.data() + std::size_t(d) * (field.bins + 1);
                const float own = float((p[hi] - p[lo]) / count);
                smoothed_.at(d, b, i) = denominator > 0.f ? total * (own + prior) / denominator : 0.f;
            }
        }
    }
}

std::uint32_t ReverbIrSynth::build(const EnergyField& field, const ReverbIrShape& shape, float* const* out) {
    for (std::uint32_t c = 0; c < kAmbisonicChannels; ++c) std::fill_n(out[c], length_, 0.f);
    if (field.bins == 0 || field.paths <= 0.0) return 0;
    smooth(field);
    const auto end = std::uint32_t(std::ceil(shape.fadeOutEnd * sampleRate_));
    const std::uint32_t samples = std::min(length_, end > shape.advance ? end - shape.advance : 0u);
    const float samplesPerBin = field.binSeconds * sampleRate_;
    // Shared by every direction and band: the window, and where each sample falls between bin centres.
    window_.resize(samples);
    position_.resize(samples);
    for (std::uint32_t n = 0; n < samples; ++n) {
        const float seconds = float(n + shape.advance) / sampleRate_;
        window_[n] = smoothstep(shape.fadeInStart, shape.fadeInEnd, seconds) * (1.f - smoothstep(shape.fadeOutStart, shape.fadeOutEnd, seconds));
        position_[n] = std::max(0.f, seconds / field.binSeconds - 0.5f);
    }
    const auto& directions = fieldDirections();
    std::vector<float>& amplitude = amplitude_;
    amplitude.resize(field.bins + 1);
    for (std::uint32_t d = 0; d < kFieldDirections; ++d) {
        std::fill_n(signal_.data(), samples, 0.f);
        const std::uint32_t shift = d * (length_ / kFieldDirections);
        for (std::uint32_t b = 0; b < kBandCount; ++b) {
            bool any = false;
            for (std::uint32_t i = 0; i < field.bins; ++i) {
                amplitude[i] = std::sqrt(std::max(smoothed_.at(d, b, i), 0.f) / samplesPerBin);
                any = any || amplitude[i] > 0.f;
            }
            if (!any) continue;
            amplitude[field.bins] = 0.f;
            for (std::uint32_t n = 0; n < samples; ++n) {
                const float position = position_[n];
                const auto i0 = std::min<std::uint32_t>(std::uint32_t(position), field.bins);
                const float frac = position - float(i0);
                const float a0 = amplitude[i0], a1 = amplitude[std::min(i0 + 1, field.bins)];
                envelope_[n] = window_[n] * (a0 + (a1 - a0) * frac);
            }
            const float* noise = noise_.data() + std::size_t(b) * length_;
            const std::uint32_t first = std::min(samples, length_ - shift);
            vDSP_vma(envelope_.data(), 1, noise + shift, 1, signal_.data(), 1, signal_.data(), 1, first);
            if (samples > first) vDSP_vma(envelope_.data() + first, 1, noise, 1, signal_.data() + first, 1, signal_.data() + first, 1, samples - first);
        }
        // AmbiX axes in the engine's world: front is -Z, left is -X, up is +Y.
        const Vec3 u = directions[d];
        const float gains[kAmbisonicChannels] = {1.f, -u.x, u.y, -u.z};
        for (std::uint32_t c = 0; c < kAmbisonicChannels; ++c) vDSP_vsma(signal_.data(), 1, &gains[c], out[c], 1, out[c], 1, samples);
    }
    return samples;
}

}  // namespace mspatial::acoustics
