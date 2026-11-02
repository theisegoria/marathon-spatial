// SPDX-License-Identifier: Apache-2.0
// Turns measured impulse responses into the engine's form. Each response splits into a minimum-phase filter and an onset delay
// (Nam, Abel and Smith 2008), because aligned minimum-phase filters interpolate without the combing that raw responses produce,
// and the interaural delay can then move continuously and separately. Off the render side only.
#include <dispatch/dispatch.h>

#include <algorithm>
#include <cmath>
#include <vector>

#include "HrtfSet.hpp"
#include "core/Dsp.hpp"
#include "core/Fft.hpp"

namespace mspatial::hrtf {

namespace {

constexpr double kSpeedOfSound = 343.0;

// Packed one-sided spectra throughout: re[0] is DC, im[0] is Nyquist.
struct Workspace {
    core::RealFft fft;
    std::vector<float> re, im, time, cepstrum;
    explicit Workspace(std::uint32_t N) : fft(N), re(N / 2), im(N / 2), time(N), cepstrum(N) {}
};

void minimumPhaseFromLog(Workspace& w, const std::vector<float>& logMagnitude, float* outRe, float* outIm) {
    const std::uint32_t N = w.fft.size();
    w.re[0] = logMagnitude[0];
    w.im[0] = logMagnitude[N / 2];
    for (std::uint32_t k = 1; k < N / 2; ++k) {
        w.re[k] = logMagnitude[k];
        w.im[k] = 0.0f;
    }
    w.fft.inverse(w.re.data(), w.im.data(), w.cepstrum.data());
    const float scale = 1.0f / static_cast<float>(N);
    w.cepstrum[0] *= scale;
    for (std::uint32_t n = 1; n < N / 2; ++n) w.cepstrum[n] *= 2.0f * scale;
    w.cepstrum[N / 2] *= scale;
    std::fill(w.cepstrum.begin() + N / 2 + 1, w.cepstrum.end(), 0.0f);
    w.fft.forward(w.cepstrum.data(), outRe, outIm);
    outRe[0] = std::exp(outRe[0] * 0.5f);
    outIm[0] = std::exp(outIm[0] * 0.5f);
    for (std::uint32_t k = 1; k < N / 2; ++k) {
        const float amplitude = std::exp(outRe[k] * 0.5f), phase = outIm[k] * 0.5f;
        outRe[k] = amplitude * std::cos(phase);
        outIm[k] = amplitude * std::sin(phase);
    }
}

std::vector<double> thirdOctaveSmooth(const std::vector<double>& power, double binHz) {
    const std::size_t count = power.size();
    std::vector<double> prefix(count + 1, 0.0), out(count);
    for (std::size_t k = 0; k < count; ++k) prefix[k + 1] = prefix[k] + power[k];
    const double half = std::pow(2.0, 1.0 / 6.0);
    for (std::size_t k = 0; k < count; ++k) {
        const double f = std::max(static_cast<double>(k), 1.0) * binHz;
        const auto lo = static_cast<std::size_t>(std::floor(f / half / binHz));
        const auto hi = std::min<std::size_t>(static_cast<std::size_t>(std::ceil(f * half / binHz)), count - 1);
        out[k] = (prefix[hi + 1] - prefix[lo]) / static_cast<double>(hi + 1 - lo);
    }
    return out;
}

}  // namespace

std::unique_ptr<HrtfSet> buildSet(const RawHrirSet& raw, const BuildOptions& options) {
    const float rate = options.sampleRate > 0 ? options.sampleRate : raw.sampleRate;
    const auto directions = static_cast<std::uint32_t>(raw.directions.size());
    if (directions == 0 || raw.taps == 0) return nullptr;

    // Bring the responses to the target rate first.
    std::vector<float> resampled;
    const float* source = raw.samples.data();
    std::uint32_t rawTaps = raw.taps;
    if (rate != raw.sampleRate) {
        rawTaps = static_cast<std::uint32_t>(core::resampledLength(raw.taps, raw.sampleRate, rate));
        resampled.resize(static_cast<std::size_t>(directions) * 2 * rawTaps);
        float* outputs = resampled.data();
        const float* inputs = raw.samples.data();
        const std::uint32_t inTaps = raw.taps, outTaps = rawTaps;
        const float inRate = raw.sampleRate;
        dispatch_apply(directions * 2u, DISPATCH_APPLY_AUTO, ^(std::size_t channel) {
            float* out = outputs + channel * outTaps;
            core::resampleSinc(inputs + channel * inTaps, inTaps, 1, inRate, rate, out, outTaps);
            for (std::uint32_t n = 0; n < outTaps; ++n) out[n] *= inRate / rate;  // Keeps the response's level: see resampleSet.
        });
        source = resampled.data();
    }

    auto set = std::make_unique<HrtfSet>();
    set->name = raw.name;
    set->licence = raw.licence;
    set->sampleRate = rate;
    set->distance = raw.distance;
    set->directions = directions;
    set->taps = options.taps ? options.taps : std::min<std::uint32_t>(512, static_cast<std::uint32_t>(std::lround(256.0 * rate / 48000.0)));
    set->dirs.reset(directions * 3ull);
    set->delays.reset(directions * 2ull);
    set->filters.reset(directions * 2ull * set->taps);
    for (std::uint32_t d = 0; d < directions; ++d) {
        set->dirs[d * 3ull] = raw.directions[d].x;
        set->dirs[d * 3ull + 1] = raw.directions[d].y;
        set->dirs[d * 3ull + 2] = raw.directions[d].z;
    }
    set->prepareLookup();

    std::uint32_t N = 2048;
    while (N < 4 * std::max(rawTaps, set->taps)) N *= 2;
    const std::uint32_t bins = N / 2 + 1;
    const double binHz = static_cast<double>(rate) / N;
    const bool extendBass = options.diffuseFieldEq;  // Loudspeakers measure little under 150 Hz, where a head is transparent.
    const auto bassBin = std::max<std::uint32_t>(2, static_cast<std::uint32_t>(150.0 / binHz));

    const auto spectrum = [&](Workspace& w, const float* ir, std::vector<float>& logMagnitude) {
        std::fill(w.time.begin(), w.time.end(), 0.0f);
        std::copy(ir, ir + rawTaps, w.time.begin());
        w.fft.forward(w.time.data(), w.re.data(), w.im.data());
        double peak = 0.0;
        logMagnitude[0] = 0.5f * std::fabs(w.re[0]);
        logMagnitude[N / 2] = 0.5f * std::fabs(w.im[0]);
        for (std::uint32_t k = 1; k < N / 2; ++k) logMagnitude[k] = 0.5f * std::hypot(w.re[k], w.im[k]);
        for (std::uint32_t k = 0; k < bins; ++k) peak = std::max<double>(peak, logMagnitude[k]);
        const float floorValue = std::max(static_cast<float>(peak) * 1e-5f, 1e-12f);
        for (std::uint32_t k = 0; k < bins; ++k) logMagnitude[k] = std::log(std::max(logMagnitude[k], floorValue));
        if (extendBass)
            for (std::uint32_t k = 0; k < bassBin; ++k) logMagnitude[k] = logMagnitude[bassBin];
    };

    // The diffuse-field average, weighted by solid angle. Its inverse, in minimum phase, is a common filter on every response.
    std::vector<float> eqRe(N / 2), eqIm(N / 2);
    if (options.diffuseFieldEq) {
        std::vector<double> power(bins, 0.0);
        Workspace w(N);
        std::vector<float> logMagnitude(bins);
        for (std::uint32_t d = 0; d < directions; ++d)
            for (std::uint32_t ear = 0; ear < 2; ++ear) {
                spectrum(w, source + (d * 2ull + ear) * rawTaps, logMagnitude);
                const double weight = set->quadrature[d] * 0.5;
                for (std::uint32_t k = 0; k < bins; ++k) {
                    const double m = std::exp(2.0 * logMagnitude[k]);
                    power[k] += weight * m;
                }
            }
        const std::vector<double> smooth = thirdOctaveSmooth(power, binHz);
        const auto top = static_cast<std::uint32_t>(std::min(18000.0, 0.45 * rate) / binHz);
        std::vector<float> logEq(bins);
        for (std::uint32_t k = 0; k < bins; ++k) {
            const double value = -0.5 * std::log(std::max(smooth[std::min(k, top)], 1e-20));
            logEq[k] = static_cast<float>(std::clamp(value, -std::log(100.0) + (-0.5 * std::log(std::max(smooth[bassBin * 2], 1e-20))),
                                                     std::log(100.0) + (-0.5 * std::log(std::max(smooth[bassBin * 2], 1e-20)))));
        }
        minimumPhaseFromLog(w, logEq, eqRe.data(), eqIm.data());
    }

    HrtfSet* target = set.get();
    const float* storedDelays = raw.delays.empty() ? nullptr : raw.delays.data();
    const float delayScale = rate / raw.sampleRate;
    const float* eqReal = eqRe.data();
    const float* eqImag = eqIm.data();
    const bool equalise = options.diffuseFieldEq;
    const auto* spectrumOf = &spectrum;
    dispatch_apply(directions, DISPATCH_APPLY_AUTO, ^(std::size_t d) {
        Workspace w(N);
        std::vector<float> logMagnitude(bins), xRe(N / 2), xIm(N / 2), mRe(N / 2), mIm(N / 2), scratchRe(N / 2), scratchIm(N / 2);
        for (std::uint32_t ear = 0; ear < 2; ++ear) {
            const float* ir = source + (d * 2 + ear) * rawTaps;
            (*spectrumOf)(w, ir, logMagnitude);
            std::copy(w.re.begin(), w.re.end(), xRe.begin());
            std::copy(w.im.begin(), w.im.end(), xIm.begin());
            minimumPhaseFromLog(w, logMagnitude, mRe.data(), mIm.data());

            // The excess-phase delay is where the response correlates best with its own minimum-phase form, over the band
            // that carries reliable timing.
            const auto correlationTop = std::min(N / 2, static_cast<std::uint32_t>(8000.0 / binHz));
            scratchRe[0] = 0.5f * xRe[0] * mRe[0];
            scratchIm[0] = 0.0f;
            for (std::uint32_t k = 1; k < N / 2; ++k) {
                const float xr = 0.5f * xRe[k], xi = 0.5f * xIm[k];
                const float keep = k < correlationTop ? 1.0f : 0.0f;
                scratchRe[k] = keep * (xr * mRe[k] + xi * mIm[k]);
                scratchIm[k] = keep * (xi * mRe[k] - xr * mIm[k]);
            }
            w.fft.inverse(scratchRe.data(), scratchIm.data(), w.time.data());
            std::uint32_t best = 0;
            float bestValue = w.time[0];
            for (std::uint32_t n = 1; n < N / 2; ++n)
                if (w.time[n] > bestValue) {
                    bestValue = w.time[n];
                    best = n;
                }
            double offset = 0.0;
            if (best > 0 && best + 1 < N) {
                const double a = w.time[best - 1], b = w.time[best], c = w.time[best + 1], denominator = a - 2.0 * b + c;
                if (std::fabs(denominator) > 1e-20) offset = 0.5 * (a - c) / denominator;
            }
            target->delays[d * 2 + ear] = static_cast<float>(best + offset) + (storedDelays ? storedDelays[d * 2 + ear] * delayScale : 0.0f);

            // The filter: minimum phase, equalised, truncated with a short fade.
            if (equalise) {
                scratchRe[0] = mRe[0] * eqReal[0];
                scratchIm[0] = mIm[0] * eqImag[0];
                for (std::uint32_t k = 1; k < N / 2; ++k) {
                    scratchRe[k] = mRe[k] * eqReal[k] - mIm[k] * eqImag[k];
                    scratchIm[k] = mRe[k] * eqImag[k] + mIm[k] * eqReal[k];
                }
            } else {
                std::copy(mRe.begin(), mRe.end(), scratchRe.begin());
                std::copy(mIm.begin(), mIm.end(), scratchIm.begin());
            }
            w.fft.inverse(scratchRe.data(), scratchIm.data(), w.time.data());
            const std::uint32_t taps = target->taps;
            float* out = target->filters.data() + (d * 2 + ear) * taps;
            const std::uint32_t fade = std::min<std::uint32_t>(32, taps / 4);
            const float scale = 1.0f / static_cast<float>(N);
            for (std::uint32_t n = 0; n < taps; ++n) {
                float value = w.time[n] * scale;
                if (n + fade >= taps) {
                    const float t = static_cast<float>(taps - n) / static_cast<float>(fade);
                    value *= std::sin(0.5f * core::kPi * t) * std::sin(0.5f * core::kPi * t);
                }
                out[n] = value;
            }
        }
    });

    float smallest = set->delays[0];
    for (std::size_t i = 0; i < directions * 2ull; ++i) smallest = std::min(smallest, set->delays[i]);
    float largest = 0.0f, widest = 0.0f;
    for (std::size_t i = 0; i < directions * 2ull; ++i) {
        set->delays[i] -= smallest;
        largest = std::max(largest, set->delays[i]);
    }
    for (std::uint32_t d = 0; d < directions; ++d) widest = std::max(widest, std::fabs(set->delays[d * 2ull] - set->delays[d * 2ull + 1]));
    set->maxDelay = largest;
    set->headRadius = options.headRadius > 0
                          ? options.headRadius
                          : std::clamp(static_cast<float>(widest / rate * kSpeedOfSound / (M_PI_2 + 1.0)), 0.06f, 0.12f);
    set->finalize();
    return set;
}

std::unique_ptr<HrtfSet> resampleSet(const HrtfSet& in, float sampleRate) {
    auto set = std::make_unique<HrtfSet>();
    set->name = in.name;
    set->licence = in.licence;
    set->sampleRate = sampleRate;
    set->headRadius = in.headRadius;
    set->distance = in.distance;
    set->directions = in.directions;
    // A minimum-phase filter has its energy at sample 0, and at a new rate that peak spreads either side of it. What falls before
    // the first sample is part of the level, so it is kept, at the cost of this constant delay on every direction.
    constexpr double kLeadSeconds = 8.0 / 48000.0;
    const auto lead = static_cast<std::uint32_t>(std::lround(kLeadSeconds * sampleRate));
    set->taps = std::min<std::uint32_t>(512, static_cast<std::uint32_t>(std::lround(in.taps * static_cast<double>(sampleRate) / in.sampleRate)) + lead);
    set->dirs.reset(in.directions * 3ull);
    std::copy(in.dirs.data(), in.dirs.data() + in.directions * 3ull, set->dirs.data());
    set->delays.reset(in.directions * 2ull);
    const float ratio = sampleRate / in.sampleRate;
    for (std::size_t i = 0; i < in.directions * 2ull; ++i) set->delays[i] = in.delays[i] * ratio;
    set->maxDelay = in.maxDelay * ratio;
    set->filters.reset(in.directions * 2ull * set->taps);
    // An impulse response holds its response in the sum of its samples, so resampled as a waveform its level moves by
    // outRate / inRate. Undo that, or a 96 kHz device hears every set 6 dB hot.
    const float gain = in.sampleRate / sampleRate;
    HrtfSet* raw = set.get();
    const HrtfSet* origin = &in;
    dispatch_apply(in.directions * 2u, DISPATCH_APPLY_AUTO, ^(std::size_t channel) {
        float* out = raw->filters.data() + channel * raw->taps;
        core::resampleSinc(origin->filters.data() + channel * origin->taps, origin->taps, 1, origin->sampleRate, sampleRate, out, raw->taps, lead);
        for (std::uint32_t n = 0; n < raw->taps; ++n) out[n] *= gain;
    });
    set->finalize();
    return set;
}

void HrtfSet::prepareLookup() {
    grid_.reset(static_cast<std::size_t>(kGrid) * kGrid * kCandidates);
    quadrature.reset(directions);
    const auto count = directions;
    const float* d = dirs.data();
    std::uint32_t* grid = grid_.data();
    dispatch_apply(kGrid, DISPATCH_APPLY_AUTO, ^(std::size_t row) {
        std::vector<std::pair<float, std::uint32_t>> ranked(count);
        for (std::uint32_t column = 0; column < kGrid; ++column) {
            // The centre of an octahedral-map cell.
            float u = (static_cast<float>(column) + 0.5f) / kGrid * 2.0f - 1.0f;
            float v = (static_cast<float>(row) + 0.5f) / kGrid * 2.0f - 1.0f;
            float w = 1.0f - std::fabs(u) - std::fabs(v);
            if (w < 0.0f) {
                const float uu = (1.0f - std::fabs(v)) * (u >= 0.0f ? 1.0f : -1.0f);
                const float vv = (1.0f - std::fabs(u)) * (v >= 0.0f ? 1.0f : -1.0f);
                u = uu;
                v = vv;
            }
            const Vec3 centre = core::normalized(Vec3{u, v, w});
            for (std::uint32_t i = 0; i < count; ++i)
                ranked[i] = {-(centre.x * d[i * 3ull] + centre.y * d[i * 3ull + 1] + centre.z * d[i * 3ull + 2]), i};
            const auto take = std::min<std::size_t>(kCandidates, count);
            std::partial_sort(ranked.begin(), ranked.begin() + static_cast<std::ptrdiff_t>(take), ranked.end());
            std::uint32_t* cell = grid + (row * kGrid + column) * kCandidates;
            for (std::uint32_t i = 0; i < kCandidates; ++i) cell[i] = ranked[std::min<std::size_t>(i, take - 1)].second;
        }
    });
    // Solid angle follows the local spacing: the squared distance to the fifth neighbour.
    double total = 0.0;
    for (std::uint32_t i = 0; i < count; ++i) {
        Neighbours n;
        const Vec3 direction{d[i * 3ull], d[i * 3ull + 1], d[i * 3ull + 2]};
        const float reach = fifthDistance(direction, n);
        quadrature[i] = reach * reach;
        total += quadrature[i];
    }
    for (std::uint32_t i = 0; i < count; ++i) quadrature[i] = static_cast<float>(quadrature[i] / total);
}

void HrtfSet::finalize() {
    if (!grid_.size()) prepareLookup();
    // Near-field shelves from the sphere model, relative to the measured distance.
    static constexpr float kRhos[kNearRhos] = {1.15f, 1.35f, 1.7f, 2.2f, 3.0f, 4.5f, 7.0f, 12.0f, 24.0f};
    const double reference = std::max(1.0, static_cast<double>(distance) / headRadius);
    constexpr double speed = kSpeedOfSound;
    constexpr int kProbes = 20;
    for (int r = 0; r < kNearRhos; ++r)
        for (int a = 0; a < kNearAngles; ++a) {
            NearFieldEntry& entry = near_[r][a];
            if (kRhos[r] >= reference) {
                entry = {0.0f, 0.0f, 1000.0f};
                continue;
            }
            const double x = std::cos(a * M_PI / (kNearAngles - 1));
            double db[kProbes];
            for (int p = 0; p < kProbes; ++p) {
                const double hz = 40.0 * std::pow(10000.0 / 40.0, static_cast<double>(p) / (kProbes - 1));
                const double ka = 2.0 * M_PI * hz * headRadius / speed;
                double reN, imN, reF, imF;
                sphereTransfer(kRhos[r], x, ka, reN, imN);
                sphereTransfer(reference, x, ka, reF, imF);
                db[p] = 10.0 * std::log10(std::max(reN * reN + imN * imN, 1e-12) / std::max(reF * reF + imF * imF, 1e-12));
            }
            const double dc = db[0], hf = db[kProbes - 1], middle = 0.5 * (dc + hf);
            double corner = 1000.0;
            for (int p = 1; p < kProbes; ++p)
                if ((db[p - 1] - middle) * (db[p] - middle) <= 0.0) {
                    const double t = (middle - db[p - 1]) / (db[p] - db[p - 1] + 1e-12);
                    corner = 40.0 * std::pow(10000.0 / 40.0, (p - 1 + t) / (kProbes - 1));
                    break;
                }
            entry = {static_cast<float>(dc), static_cast<float>(hf), static_cast<float>(corner)};
        }
}

float HrtfSet::fifthDistance(Vec3 direction, Neighbours& out) const noexcept {
    const float ax = std::fabs(direction.x), ay = std::fabs(direction.y), az = std::fabs(direction.z);
    const float s = 1.0f / std::max(ax + ay + az, 1e-9f);
    float u = direction.x * s, v = direction.y * s;
    if (direction.z < 0.0f) {
        const float uu = (1.0f - std::fabs(v)) * (u >= 0.0f ? 1.0f : -1.0f);
        const float vv = (1.0f - std::fabs(u)) * (v >= 0.0f ? 1.0f : -1.0f);
        u = uu;
        v = vv;
    }
    const auto column = std::min<std::uint32_t>(kGrid - 1, static_cast<std::uint32_t>((u * 0.5f + 0.5f) * kGrid));
    const auto row = std::min<std::uint32_t>(kGrid - 1, static_cast<std::uint32_t>((v * 0.5f + 0.5f) * kGrid));
    const std::uint32_t* cell = grid_.data() + (static_cast<std::size_t>(row) * kGrid + column) * kCandidates;
    // Insertion into the five nearest. Candidates repeat only in sets smaller than kCandidates.
    float bestDot[5] = {-2.0f, -2.0f, -2.0f, -2.0f, -2.0f};
    std::uint32_t bestIndex[5] = {};
    const float* d = dirs.data();
    std::uint32_t previous = 0xffffffffu;
    for (std::uint32_t i = 0; i < kCandidates; ++i) {
        const std::uint32_t index = cell[i];
        if (index == previous) continue;
        previous = index;
        const float value = direction.x * d[index * 3ull] + direction.y * d[index * 3ull + 1] + direction.z * d[index * 3ull + 2];
        if (value <= bestDot[4]) continue;
        int slot = 4;
        while (slot > 0 && bestDot[slot - 1] < value) {
            bestDot[slot] = bestDot[slot - 1];
            bestIndex[slot] = bestIndex[slot - 1];
            --slot;
        }
        bestDot[slot] = value;
        bestIndex[slot] = index;
    }
    for (int i = 0; i < 4; ++i) out.index[i] = bestIndex[i];
    for (int i = 0; i < 4; ++i) out.weight[i] = bestDot[i];
    return std::sqrt(std::max(0.0f, 2.0f - 2.0f * std::max(bestDot[4], -1.0f)));
}

void HrtfSet::query(Vec3 direction, Neighbours& out) const noexcept {
    const float radius = fifthDistance(direction, out);
    float total = 0.0f;
    for (int i = 0; i < 4; ++i) {
        const float cosine = out.weight[i];
        if (cosine < -1.5f) {  // Fewer than four measurements in the set.
            out.weight[i] = 0.0f;
            continue;
        }
        const float distance = std::sqrt(std::max(0.0f, 2.0f - 2.0f * cosine));
        if (distance < 1e-5f) {
            const std::uint32_t only = out.index[i];
            out.index[0] = out.index[1] = out.index[2] = out.index[3] = only;
            out.weight[0] = 1.0f;
            out.weight[1] = out.weight[2] = out.weight[3] = 0.0f;
            return;
        }
        const float share = std::max(radius - distance, 0.0f) / (radius * distance + 1e-12f);
        out.weight[i] = share * share;
        total += out.weight[i];
    }
    if (total <= 0.0f) {
        out.weight[0] = 1.0f;
        out.weight[1] = out.weight[2] = out.weight[3] = 0.0f;
        return;
    }
    for (int i = 0; i < 4; ++i) out.weight[i] /= total;
}

NearFieldEntry HrtfSet::nearField(float rho, float angle) const noexcept {
    static constexpr float kRhos[kNearRhos] = {1.15f, 1.35f, 1.7f, 2.2f, 3.0f, 4.5f, 7.0f, 12.0f, 24.0f};
    const float clamped = std::clamp(rho, kRhos[0], kRhos[kNearRhos - 1]);
    int r = 0;
    while (r < kNearRhos - 2 && clamped > kRhos[r + 1]) ++r;
    const float tr = (std::log(clamped) - std::log(kRhos[r])) / (std::log(kRhos[r + 1]) - std::log(kRhos[r]));
    const float position = std::clamp(angle / core::kPi, 0.0f, 1.0f) * (kNearAngles - 1);
    const int a = std::min(static_cast<int>(position), kNearAngles - 2);
    const float ta = position - static_cast<float>(a);
    const auto mix = [&](float NearFieldEntry::*field) {
        const float low = near_[r][a].*field * (1 - ta) + near_[r][a + 1].*field * ta;
        const float high = near_[r + 1][a].*field * (1 - ta) + near_[r + 1][a + 1].*field * ta;
        return low * (1 - tr) + high * tr;
    };
    return {mix(&NearFieldEntry::dcDb), mix(&NearFieldEntry::hfDb), mix(&NearFieldEntry::cornerHz)};
}

}  // namespace mspatial::hrtf
