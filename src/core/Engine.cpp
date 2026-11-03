// SPDX-License-Identifier: Apache-2.0
// The render graph: one block at a time, no locks, no allocation, no clock in the signal path.
#include "Engine.hpp"

#include <mach/mach_time.h>

#include <algorithm>
#include <cmath>

using namespace mspatial::core;
using mspatial::hrtf::HrtfSet;
using mspatial::hrtf::Neighbours;
using mspatial::hrtf::NearFieldEntry;

namespace mspatial::core {
namespace {

constexpr float kBaseDelay = 4.0f;               // The interpolator looks four samples ahead.
constexpr float kFilterThreshold = 0.99999f;     // Resynthesise past about 0.25 degrees.
constexpr std::uint32_t kReflectionOrder = 3;
constexpr std::uint32_t kRankInterval = 8;
constexpr float kVoiceHoldMargin = 1.5f;         // About 3.5 dB: a voice in use is only given up to something clearly louder.
constexpr float kMinPitch = 1.0f / 32.0f;        // 8 kHz content on a 192 kHz device, with room for Doppler below it.
constexpr float kMinRate = 1.0f / 64.0f;

inline float dbToGain(float db) noexcept { return std::pow(10.0f, db / 20.0f); }

inline bool finiteBlock(const float* x, std::uint32_t n) noexcept {
    float energy = 0.0f;
    vDSP_svesq(x, 1, &energy, n);  // NaN and infinity survive a sum of squares; a maximum can skip over a NaN.
    return std::isfinite(energy);
}

// Blocks the ambisonic mix keeps decoding after its last input: the binaural decoder's filter length, or a block or two of a speaker matrix.
inline std::uint32_t mixHoldBlocks(const ms_engine& e) noexcept { return e.speakers ? 2 : e.decoder ? e.decoder->partitions + 2 : 0; }

inline Quat unpackQuat(std::uint64_t v) noexcept {
    const auto unpack = [](std::uint64_t bits) { return static_cast<float>(static_cast<std::int16_t>(bits & 0xffffu)) / 32767.0f; };
    return normalized(Quat{unpack(v), unpack(v >> 16), unpack(v >> 32), unpack(v >> 48)});
}

float distanceGain(const ms_source_params& p, float d) noexcept {
    const float ref = std::max(p.ref_distance, 1e-3f);
    switch (p.distance_model) {
        case MS_DISTANCE_INVERSE: {
            const float c = std::clamp(d, ref, std::max(p.max_distance, ref));
            return ref / (ref + p.rolloff * (c - ref));
        }
        case MS_DISTANCE_LINEAR: {
            const float top = std::max(p.max_distance, ref + 1e-3f), c = std::clamp(d, ref, top);
            return std::clamp(1.0f - p.rolloff * (c - ref) / (top - ref), 0.0f, 1.0f);
        }
        case MS_DISTANCE_EXPONENTIAL: {
            const float c = std::clamp(d, ref, std::max(p.max_distance, ref));
            return std::pow(c / ref, -p.rolloff);
        }
        default: return 1.0f;
    }
}

// Ramps the block in place from a to b.
inline void applyRamp(float* x, std::uint32_t n, float a, float b) noexcept {
    if (a == b) {
        if (a != 1.0f) vDSP_vsmul(x, 1, &a, x, 1, n);
        return;
    }
    const float step = (b - a) / static_cast<float>(n);
    for (std::uint32_t i = 0; i < n; ++i) x[i] *= a + step * static_cast<float>(i + 1);
}

inline void cmacPacked(float* __restrict ar, float* __restrict ai, const float* xr, const float* xi, const float* hr, const float* hi, std::uint32_t bins) noexcept {
    ar[0] += xr[0] * hr[0];
    ai[0] += xi[0] * hi[0];
#pragma clang loop vectorize(enable)
    for (std::uint32_t k = 1; k < bins; ++k) {
        ar[k] += xr[k] * hr[k] - xi[k] * hi[k];
        ai[k] += xr[k] * hi[k] + xi[k] * hr[k];
    }
}

inline void interpolateNeighbours(const HrtfSet& set, const Neighbours& n, std::uint32_t ear, float* out) noexcept {
    const std::uint32_t taps = set.taps;
    std::memset(out, 0, taps * sizeof(float));
    for (int i = 0; i < 4; ++i) {
        if (n.weight[i] == 0.0f) continue;
        const float* h = set.filters.data() + (n.index[i] * 2ull + ear) * taps;
        vDSP_vsma(h, 1, &n.weight[i], out, 1, out, 1, taps);
    }
}

inline float interpolateDelay(const HrtfSet& set, const Neighbours& n, std::uint32_t ear) noexcept {
    float sum = 0.0f;
    for (int i = 0; i < 4; ++i) sum += n.weight[i] * set.delays[n.index[i] * 2ull + ear];
    return sum;
}

}  // namespace
}  // namespace mspatial::core

namespace {

void clearBiquads(RenderSource& s) {
    for (Biquad& q : s.stages) q = Biquad{};
    for (float& v : s.stageValue) v = NAN;
}

void resetSource(RenderSource& s) {
    s.state = RenderSource::State::Idle;
    s.buffer = nullptr;
    s.channel = 0;
    s.loop = false;
    s.callback = nullptr;
    s.user = nullptr;
    s.position = 0.0;
    s.seekTarget = 0.0;
    s.scheduled = false;
    s.primed = false;
    s.gain = 0.0f;
    s.dryGain = 0.0f;
    s.rate = 1.0f;
    s.pan[0] = s.pan[1] = 0.0f;
    s.pullFill = 0;
    s.pullFraction = 0.0f;
    s.pullEnded = false;
    s.contentEnded = false;
    s.tailBlocks = 0;
    s.directMix = 0.0f;
    s.wantDirect = false;
    s.score = 0.0f;
    s.encodeWeight = 0.0f;
    s.peak = 0.0f;
    std::memset(s.encode, 0, sizeof s.encode);
    s.haveSpeakerGain = false;
    clearBiquads(s);
}

void releaseVoice(ms_engine& e, RenderSource& s) {
    if (!s.voice) return;
    e.freeVoices[e.freeVoiceCount++] = static_cast<std::uint32_t>(s.voice - e.voices.data());
    s.voice = nullptr;
    s.voiceDrain = 0;
}

// How long a voice's state keeps sounding once nothing more goes in.
std::uint32_t voiceDrainBlocks(const ms_engine& e, const HrtfVoice& v) {
    const std::uint32_t B = e.B;
    std::uint32_t frames = e.set->taps + static_cast<std::uint32_t>(std::ceil(e.set->maxDelay)) + B;
    float longest = 0.0f;
    for (const HrtfVoice::Tap& t : v.taps)
        if (t.active) longest = std::max(longest, t.delay);
    frames += static_cast<std::uint32_t>(std::ceil(longest));
    return (frames + B - 1) / B + 2;
}

void resetVoiceState(ms_engine& e, HrtfVoice& v) {
    const std::uint32_t B = e.B;
    std::memset(v.input, 0, 2ull * B * sizeof(float));
    std::memset(v.fdlRe, 0, static_cast<std::size_t>(e.maxPartitions) * B * sizeof(float));
    std::memset(v.fdlIm, 0, static_cast<std::size_t>(e.maxPartitions) * B * sizeof(float));
    std::memset(v.delay[0], 0, 2ull * kDelayFrames * sizeof(float));
    std::memset(v.delay[1], 0, 2ull * kDelayFrames * sizeof(float));
    for (float* band : v.bands) std::memset(band, 0, kReflectionFrames * sizeof(float));
    v.fdlHead = 0;
    v.delayWrite = 0;
    v.bandWrite = 0;
    v.delayNow[0] = v.delayNow[1] = 0.0f;
    v.haveFilter = false;
    v.crossfade = false;
    v.lastSet = nullptr;
    v.near[0] = v.near[1] = Biquad{};
    v.nearGain[0] = v.nearGain[1] = 1.0f;
    v.splitLow.clear();
    v.splitHigh.clear();
    for (HrtfVoice::Tap& t : v.taps) t = HrtfVoice::Tap{};
}

void acquireVoice(ms_engine& e, RenderSource& s) {
    if (s.voice || e.freeVoiceCount == 0) return;
    s.voice = &e.voices[e.freeVoices[--e.freeVoiceCount]];
    resetVoiceState(e, *s.voice);
    s.directMix = 0.0f;
}

void processCommands(ms_engine& e) {
    Command c;
    while (e.commands.pop(c)) {
        switch (c.type) {
            case CommandType::SourceReset: {
                RenderSource& s = e.sources[c.index];
                releaseVoice(e, s);
                resetSource(s);
                e.published[c.index].playing.store(0, std::memory_order_relaxed);
                break;
            }
            case CommandType::SourceBuffer: {
                RenderSource& s = e.sources[c.index];
                s.buffer = static_cast<BufferObject*>(c.pointer);
                s.channel = c.a;
                s.loop = c.b != 0;
                s.callback = nullptr;
                break;
            }
            case CommandType::SourceCallback: {
                RenderSource& s = e.sources[c.index];
                s.callback = reinterpret_cast<ms_source_callback>(c.pointer);
                s.user = c.pointer2;
                s.buffer = nullptr;
                break;
            }
            case CommandType::SourcePlay: {
                RenderSource& s = e.sources[c.index];
                if (!s.buffer && !s.callback) break;
                s.state = RenderSource::State::Playing;
                s.scheduled = true;
                s.startFrame = c.frame;
                s.position = 0.0;
                s.pullFill = 0;
                s.pullFraction = 0.0f;
                s.pullEnded = false;
                s.contentEnded = false;
                s.primed = false;
                s.gain = 0.0f;
                s.dryGain = 0.0f;
                clearBiquads(s);
                break;
            }
            case CommandType::SourceStop: {
                RenderSource& s = e.sources[c.index];
                if (s.state == RenderSource::State::Playing) s.state = RenderSource::State::Stopping;
                break;
            }
            case CommandType::SourceSeek: {
                RenderSource& s = e.sources[c.index];
                if (s.scheduled || s.state == RenderSource::State::Idle) s.seekTarget = c.position;
                else s.position = c.position;
                break;
            }
            case CommandType::BedReset: {
                RenderBed& b = e.beds[c.index];
                b.state = RenderBed::State::Idle;
                b.buffer = nullptr;
                b.callback = nullptr;
                b.scheduled = false;
                b.haveRotation = false;
                b.haveGain = false;
                break;
            }
            case CommandType::BedBuffer: {
                RenderBed& b = e.beds[c.index];
                b.buffer = static_cast<BufferObject*>(c.pointer);
                b.callback = nullptr;
                b.loop = c.b != 0;
                b.order = c.a;
                break;
            }
            case CommandType::BedCallback: {
                RenderBed& b = e.beds[c.index];
                b.callback = reinterpret_cast<ms_bed_callback>(c.pointer);
                b.user = c.pointer2;
                b.buffer = nullptr;
                b.order = c.a;
                break;
            }
            case CommandType::BedPlay: {
                RenderBed& b = e.beds[c.index];
                if (!b.buffer && !b.callback) break;
                b.state = RenderBed::State::Playing;
                b.scheduled = true;
                b.startFrame = c.frame;
                b.position = 0;
                b.haveGain = false;
                break;
            }
            case CommandType::BedStop: {
                RenderBed& b = e.beds[c.index];
                if (b.state == RenderBed::State::Playing) b.state = RenderBed::State::Stopping;
                break;
            }
            case CommandType::BusProcessor: {
                RenderBus& b = e.buses[c.index];
                b.processor = reinterpret_cast<ms_bus_processor>(c.pointer);
                b.user = c.pointer2;
                break;
            }
            case CommandType::BlockCallback:
                e.blockCallback = reinterpret_cast<ms_block_callback>(c.pointer);
                e.blockUser = c.pointer2;
                break;
            case CommandType::BufferRetire: {
                auto* buffer = static_cast<BufferObject*>(c.pointer);
                for (std::uint32_t i = 0; i < e.desc.max_sources; ++i)
                    if (e.sources[i].buffer == buffer) {
                        e.sources[i].buffer = nullptr;
                        if (e.sources[i].state == RenderSource::State::Playing) e.sources[i].state = RenderSource::State::Stopping;
                        e.sources[i].contentEnded = true;
                    }
                for (std::uint32_t i = 0; i < e.desc.max_beds; ++i)
                    if (e.beds[i].buffer == buffer) {
                        e.beds[i].buffer = nullptr;
                        e.beds[i].state = RenderBed::State::Idle;
                    }
                (void)e.retired.push({Retired::Kind::Buffer, buffer, nullptr});
                break;
            }
            case CommandType::SetHrtf: {
                // The old set stays alive through this block's crossfade, and goes back once the block ends.
                if (e.retiringSet) (void)e.retired.push({Retired::Kind::Hrtf, e.retiringSet, nullptr});
                e.retiringSet = const_cast<HrtfSet*>(e.set);
                e.previousDecoder = e.decoder;
                e.set = static_cast<const HrtfSet*>(c.pointer);
                e.decoder = static_cast<const AmbisonicDecoder*>(c.pointer2);
                break;
            }
        }
    }
}

// Reads a mono block at a varying rate from a buffer. position is in frames of content.
void readBuffer(ms_engine& e, RenderSource& s, float* out, float rate0, float rate1) {
    const std::uint32_t B = e.B;
    const BufferObject* buffer = s.buffer;
    if (!buffer || buffer->frames == 0) {
        std::memset(out, 0, B * sizeof(float));
        s.contentEnded = true;
        return;
    }
    const std::uint32_t stride = buffer->channels;
    const float* data = buffer->data.data() + std::min(s.channel, stride - 1);
    const auto frames = static_cast<std::int64_t>(buffer->frames);
    double p = s.position;
    const double step = static_cast<double>(rate1 - rate0) / B;
    double rate = rate0;
    const auto at = [&](std::int64_t i) -> float {
        if (s.loop) {
            i %= frames;
            if (i < 0) i += frames;
            return data[static_cast<std::size_t>(i) * stride];
        }
        return i >= 0 && i < frames ? data[static_cast<std::size_t>(i) * stride] : 0.0f;
    };
    for (std::uint32_t t = 0; t < B; ++t) {
        rate += step;
        if (p < -4.0) {
            out[t] = 0.0f;
        } else {
            const double whole = std::floor(p);
            const auto i = static_cast<std::int64_t>(whole);
            const auto fraction = static_cast<float>(p - whole);
            if (fraction == 0.0f) {
                out[t] = at(i);
            } else {
                float window[8];
                if (i >= 3 && i + 4 < frames && stride == 1) {
                    out[t] = sincSample(data + i, fraction);
                } else {
                    for (int k = 0; k < 8; ++k) window[k] = at(i + k - 3);
                    out[t] = sincSample(window + 3, fraction);
                }
            }
            e.opReal += fraction == 0.0f ? 0 : 16;
        }
        p += rate;
        if (s.loop && p >= static_cast<double>(frames)) p -= static_cast<double>(frames);
    }
    s.position = p;
    if (!s.loop && p >= static_cast<double>(frames) + 4.0) s.contentEnded = true;
}

// Reads a mono block from a callback, through the same interpolator.
void readCallback(ms_engine& e, RenderSource& s, float* out, float rate0, float rate1) {
    const std::uint32_t B = e.B;
    constexpr std::uint32_t kHistory = 8;  // Frames behind the read head, for the interpolator.
    std::uint32_t t = 0;
    if (s.position < 0.0) {  // Scheduled inside this block, or held back by propagation delay.
        const auto lead = std::min<std::uint32_t>(B, static_cast<std::uint32_t>(std::ceil(-s.position)));
        std::memset(out, 0, lead * sizeof(float));
        s.position += lead;
        t = lead;
        if (t == B) return;
    }
    const std::uint32_t remaining = B - t;
    const std::uint32_t capacity = 5 * B + 32;
    if (rate0 == 1.0f && rate1 == 1.0f && s.pullFraction == 0.0f && s.pullFill == 0) {
        const std::uint32_t got = s.pullEnded ? 0 : s.callback(s.user, out + t, remaining);
        if (got < remaining) {
            std::memset(out + t + got, 0, (remaining - got) * sizeof(float));
            s.pullEnded = true;
            s.contentEnded = true;
        }
        s.position += got;
        return;
    }
    if (s.pullFill == 0) {  // Interpolation begins: silence behind the first sample.
        std::memset(s.pull, 0, kHistory * sizeof(float));
        s.pullFill = kHistory;
    }
    const double step = static_cast<double>(rate1 - rate0) / B;
    double total = 0.0, rate = rate0;
    for (std::uint32_t i = 0; i < remaining; ++i) {
        rate += step;
        total += rate;
    }
    // Frames from the read head to the end of this block, plus the interpolator's look-ahead.
    const auto need = static_cast<std::uint32_t>(std::ceil(s.pullFraction + total)) + 5;
    if (s.pullFill < kHistory + need) {
        const std::uint32_t want = std::min(kHistory + need - s.pullFill, capacity - s.pullFill);
        const std::uint32_t got = s.pullEnded ? 0 : s.callback(s.user, s.pull + s.pullFill, want);
        if (got < want) {
            std::memset(s.pull + s.pullFill + got, 0, (want - got) * sizeof(float));
            s.pullEnded = true;
        }
        s.pullFill += want;
    }
    double p = s.pullFraction;
    rate = rate0;
    for (; t < B; ++t) {
        rate += step;
        const double whole = std::floor(p);
        const std::uint32_t i = kHistory + static_cast<std::uint32_t>(whole);
        out[t] = i + 4 < s.pullFill ? sincSample(s.pull + i, static_cast<float>(p - whole)) : 0.0f;
        p += rate;
    }
    e.opReal += 16ull * remaining;
    s.position += p - s.pullFraction;
    const auto advance = std::min<std::uint32_t>(static_cast<std::uint32_t>(std::floor(p)), s.pullFill - kHistory);
    std::memmove(s.pull, s.pull + advance, (s.pullFill - advance) * sizeof(float));
    s.pullFill -= advance;
    s.pullFraction = static_cast<float>(p - std::floor(p));
    if (s.pullEnded) s.contentEnded = true;  // The read head has passed everything the source gave.
}

// Updates one filter stage if its parameter moved, in steps inside the block so a sweep stays smooth.
template <class Design>
inline void updateStage(RenderSource& s, int stage, float value, float tolerance, Design&& design) {
    if (s.stageValue[stage] == s.stageValue[stage] && std::fabs(value - s.stageValue[stage]) < tolerance) return;
    s.stageValue[stage] = value;
    s.stages[stage].setCoefficients(design());
}

void processStages(ms_engine& e, RenderSource& s, float* x) {
    for (Biquad& q : s.stages) {
        if (q.identityCoefficients() && q.silentState()) continue;
        q.process(x, e.B);
        e.opReal += 5ull * e.B;
    }
}

void synthesiseFilter(ms_engine& e, HrtfVoice& v, const Neighbours& left, const Neighbours& right) {
    const std::uint32_t B = e.B;
    const HrtfSet& set = *e.set;
    const std::uint32_t partitions = (set.taps + B - 1) / B;
    const std::uint32_t slot = 1 - v.current;
    float* time = e.scratch.data() + 8ull * B;      // [512]
    float* padded = time + 512;                     // [2 B]
    for (std::uint32_t ear = 0; ear < 2; ++ear) {
        interpolateNeighbours(set, ear == 0 ? left : right, ear, time);
        for (std::uint32_t p = 0; p < partitions; ++p) {
            const std::uint32_t begin = p * B, count = std::min(B, set.taps - begin);
            std::memcpy(padded, time + begin, count * sizeof(float));
            std::memset(padded + count, 0, (2 * B - count) * sizeof(float));
            const std::size_t at = (static_cast<std::size_t>(slot) * 2 + ear) * e.maxPartitions * B + static_cast<std::size_t>(p) * B;
            e.fft.forward(padded, v.filterRe + at, v.filterIm + at);
        }
    }
    e.opForward += 2ull * partitions;
    e.opReal += 8ull * set.taps;
    ++e.opUpdates;
    v.partitions[slot] = partitions;
    v.current = slot;
}

void convolve(ms_engine& e, const HrtfVoice& v, std::uint32_t slot, std::uint32_t ear, float* out) {
    const std::uint32_t B = e.B, P = v.partitions[slot];
    float* accRe = e.scratch.data() + 2ull * B;
    float* accIm = accRe + B;
    std::memset(accRe, 0, 2ull * B * sizeof(float));
    const float* hr = v.filterRe + (static_cast<std::size_t>(slot) * 2 + ear) * e.maxPartitions * B;
    const float* hi = v.filterIm + (static_cast<std::size_t>(slot) * 2 + ear) * e.maxPartitions * B;
    for (std::uint32_t p = 0; p < P; ++p) {
        const std::uint32_t at = (v.fdlHead + e.maxPartitions - p) % e.maxPartitions;
        cmacPacked(accRe, accIm, v.fdlRe + static_cast<std::size_t>(at) * B, v.fdlIm + static_cast<std::size_t>(at) * B,
                   hr + static_cast<std::size_t>(p) * B, hi + static_cast<std::size_t>(p) * B, B);
    }
    e.opSpectrum += static_cast<std::uint64_t>(P) * B;
    e.fft.inverse(accRe, accIm, out);
    ++e.opInverse;
}

// One ear's delay line: writes the block, reads with a delay moving linearly across the block, adds to the bus.
void delayEar(ms_engine& e, HrtfVoice& v, std::uint32_t ear, const float* in, float target, float* bus) {
    const std::uint32_t B = e.B;
    float* line = v.delay[ear];
    const std::uint32_t write = v.delayWrite;
    for (std::uint32_t t = 0; t < B; ++t) {
        const std::uint32_t at = (write + t) & (kDelayFrames - 1);
        line[at] = in[t];
        line[at + kDelayFrames] = in[t];
    }
    const float from = v.delayNow[ear];
    const float limit = static_cast<float>(kDelayFrames - B - 16) - kBaseDelay;
    const float to = std::min(target, limit);
    const float step = (to - from) / static_cast<float>(B);
    for (std::uint32_t t = 0; t < B; ++t) {
        const float d = from + step * static_cast<float>(t + 1) + kBaseDelay;
        const float whole = std::ceil(d);
        const float fraction = whole - d;  // Position = write + t - whole + fraction
        const auto index = (write + t + 2 * kDelayFrames - static_cast<std::uint32_t>(whole)) & (kDelayFrames - 1);
        bus[t] += fraction == 0.0f ? line[index] : sincSample(line + index + (index < 4 ? kDelayFrames : 0), fraction);
    }
    e.opReal += 16ull * B;
    v.delayNow[ear] = to;
}

constexpr float kMaxTapLevel = 1e6f;  // Far above any mix, and 64 taps at this level cannot overflow it.

void renderReflections(ms_engine& e, HrtfVoice& v, const SourceSnapshot& snap, const float* dry, float product) {
    const std::uint32_t B = e.B;
    float* low = e.scratch.data();
    float* high = low + B;
    bool anyActive = snap.reflectionCount > 0;
    for (const HrtfVoice::Tap& t : v.taps) anyActive = anyActive || t.active;
    if (!anyActive) return;
    std::memcpy(low, dry, B * sizeof(float));
    std::memcpy(high, dry, B * sizeof(float));
    v.splitLow.process(low, B);
    v.splitHigh.process(high, B);
    for (std::uint32_t t = 0; t < B; ++t) {
        const std::uint32_t at = (v.bandWrite + t) & (kReflectionFrames - 1);
        v.bands[0][at] = low[t];
        v.bands[2][at] = high[t];
        v.bands[1][at] = dry[t] - low[t] - high[t];
    }
    e.opReal += 12ull * B;
    const std::uint32_t order = std::min(e.desc.ambisonic_order, kReflectionOrder), channels = (order + 1) * (order + 1);
    if (channels <= 1 && order == 0) {
        v.bandWrite = (v.bandWrite + B) & (kReflectionFrames - 1);
        return;
    }
    float* mixed = e.scratch.data();
    const float invB = 1.0f / static_cast<float>(B);
    for (std::uint32_t i = 0; i < MS_MAX_REFLECTION_TAPS; ++i) {
        HrtfVoice::Tap& tap = v.taps[i];
        const bool used = i < snap.reflectionCount;
        if (!used && !tap.active) continue;
        float gainTarget[3] = {0, 0, 0};
        float encodeTarget[16] = {};
        float delayTarget = tap.delay;
        if (used) {
            const ms_reflection_tap& input = snap.reflections[i];
            for (int b = 0; b < 3; ++b) gainTarget[b] = input.gain[b] * product;  // Into the ambisonic mix, outside the bus's own signal.
            delayTarget = std::clamp(input.delay_seconds * e.rate, 1.0f, static_cast<float>(kReflectionFrames - B - 8));
            const Vec3 world = normalized(Vec3(input.direction));
            const Vec3 rel{dot(world, e.right), dot(world, e.up), dot(world, e.back)};
            encodeDirection(headToAmbisonic(rotate(e.headFrom, rel)), order, encodeTarget);
        } else {
            std::memcpy(encodeTarget, tap.encode, sizeof encodeTarget);
        }
        if (!tap.active) {  // A new path starts at its own delay, silently.
            tap.delay = delayTarget;
            std::memcpy(tap.encode, encodeTarget, sizeof encodeTarget);
            tap.gain[0] = tap.gain[1] = tap.gain[2] = 0.0f;
            tap.active = true;
        }
        const float delayStep = (delayTarget - tap.delay) * invB;
        const float g0 = tap.gain[0], g1 = tap.gain[1], g2 = tap.gain[2];
        const float s0 = (gainTarget[0] - g0) * invB, s1 = (gainTarget[1] - g1) * invB, s2 = (gainTarget[2] - g2) * invB;
        for (std::uint32_t t = 0; t < B; ++t) {
            const float ft = static_cast<float>(t);
            const float d = tap.delay + delayStep * ft;
            const float whole = std::ceil(d), fraction = whole - d;
            const std::uint32_t a = (v.bandWrite + t + kReflectionFrames - static_cast<std::uint32_t>(whole)) & (kReflectionFrames - 1);
            const std::uint32_t b = (a + 1) & (kReflectionFrames - 1);
            const float lo = v.bands[0][a] + fraction * (v.bands[0][b] - v.bands[0][a]);
            const float mi = v.bands[1][a] + fraction * (v.bands[1][b] - v.bands[1][a]);
            const float hi = v.bands[2][a] + fraction * (v.bands[2][b] - v.bands[2][a]);
            mixed[t] = (g0 + s0 * ft) * lo + (g1 + s1 * ft) * mi + (g2 + s2 * ft) * hi;
        }
        // The ambisonic mix is shared, so a tap is checked before it. A limit, not only finiteness: finite taps can still overflow the sum.
        float tapLevel = 0.0f;
        vDSP_maxmgv(mixed, 1, &tapLevel, B);
        if (!finiteBlock(mixed, B) || tapLevel > kMaxTapLevel) {
            std::memset(mixed, 0, B * sizeof(float));
            gainTarget[0] = gainTarget[1] = gainTarget[2] = 0.0f;
            e.stats.faults.fetch_add(1, std::memory_order_relaxed);
        }
        for (std::uint32_t c = 0; c < channels; ++c) {
            const float from = tap.encode[c], step = (encodeTarget[c] - from) * invB;
            float* mix = e.shMix.data() + static_cast<std::size_t>(c) * B;
#pragma clang loop vectorize(enable)
            for (std::uint32_t t = 0; t < B; ++t) mix[t] += mixed[t] * (from + step * static_cast<float>(t));
        }
        e.opReal += (9ull + channels) * B;
        tap.delay = delayTarget;
        tap.gain[0] = gainTarget[0];
        tap.gain[1] = gainTarget[1];
        tap.gain[2] = gainTarget[2];
        std::memcpy(tap.encode, encodeTarget, sizeof encodeTarget);
        if (!used) tap.active = false;
    }
    e.shActiveBlocks = mixHoldBlocks(e);
    v.bandWrite = (v.bandWrite + B) & (kReflectionFrames - 1);
}

void renderDirect(ms_engine& e, RenderSource& s, const SourceSnapshot& snap, const float* mono, float mixFrom, float mixTo, float distance, Vec3 direction) {
    const std::uint32_t B = e.B;
    HrtfVoice& v = *s.voice;
    const ms_source_params& p = snap.params;
    const HrtfSet& set = *e.set;

    // Overlap-save: the block joins the second half of the window.
    std::memcpy(v.input, v.input + B, B * sizeof(float));
    std::memcpy(v.input + B, mono, B * sizeof(float));
    applyRamp(v.input + B, B, mixFrom, mixTo);
    v.fdlHead = (v.fdlHead + 1) % e.maxPartitions;
    e.fft.forward(v.input, v.fdlRe + static_cast<std::size_t>(v.fdlHead) * B, v.fdlIm + static_cast<std::size_t>(v.fdlHead) * B);
    ++e.opForward;

    // Parallax: close to the head each ear sees its own direction. The two merge smoothly further out.
    const float radius = set.headRadius;
    Vec3 dirL = direction, dirR = direction;
    const bool near = p.near_field > 0.0f && distance < set.distance;
    if (near) {
        // The set already holds the parallax of its measured distance, so the offset fades out as the source reaches it.
        const Vec3 position = direction * std::max(distance, radius * 1.1f);
        const float offset = radius * p.near_field * (1.0f - distance / set.distance);
        dirL = normalized(position - Vec3{-offset, 0.0f, 0.0f});
        dirR = normalized(position - Vec3{offset, 0.0f, 0.0f});
    }
    const bool firstBlock = !v.haveFilter;
    if (!v.haveFilter || v.lastSet != &set || dot(dirL, v.lastDirection[0]) < kFilterThreshold || dot(dirR, v.lastDirection[1]) < kFilterThreshold) {
        Neighbours nL, nR;
        set.query(dirL, nL);
        if (near) set.query(dirR, nR);
        else nR = nL;
        const bool first = !v.haveFilter;
        synthesiseFilter(e, v, nL, nR);
        v.crossfade = !first;
        v.haveFilter = true;
        v.lastSet = &set;
        v.lastDirection[0] = dirL;
        v.lastDirection[1] = dirR;
        v.delayNow[0] = first ? interpolateDelay(set, nL, 0) : v.delayNow[0];
        v.delayNow[1] = first ? interpolateDelay(set, nR, 1) : v.delayNow[1];
        v.delayTarget[0] = interpolateDelay(set, nL, 0);
        v.delayTarget[1] = interpolateDelay(set, nR, 1);
    } else {
        v.crossfade = false;
    }

    const float scale = e.inverseScale;
    float* newer = e.scratch.data() + 4ull * B;  // [2 B]
    float* older = e.scratch.data() + 6ull * B;  // [2 B]
    float* earOut = e.scratch.data();            // [B]
    RenderBus& bus = e.buses[std::min<std::uint32_t>(p.bus ? (p.bus & 0xffffu) - 1 : 0, e.desc.max_buses - 1)];
    for (std::uint32_t ear = 0; ear < 2; ++ear) {
        convolve(e, v, v.current, ear, newer);
        if (v.crossfade) {
            convolve(e, v, 1 - v.current, ear, older);
            const float* ramp = e.crossfadeRamp.data();
#pragma clang loop vectorize(enable)
            for (std::uint32_t t = 0; t < B; ++t) earOut[t] = scale * (newer[B + t] * ramp[t] + older[B + t] * (1.0f - ramp[t]));
            e.opReal += 2ull * B;
        } else {
            vDSP_vsmul(newer + B, 1, &scale, earOut, 1, B);
            e.opReal += B;
        }
        // The near-field shelf from the sphere model; gain on the high band, a shelf for what rises below it.
        float gainTo = 1.0f;
        bool shelf = false;
        if (near) {
            const float angle = std::acos(std::clamp(direction.x * (ear == 0 ? -1.0f : 1.0f), -1.0f, 1.0f));
            const NearFieldEntry entry = set.nearField(distance / radius, angle);
            gainTo = dbToGain(entry.hfDb * p.near_field);
            v.near[ear].setCoefficients(designFirstOrderLowShelf(entry.cornerHz, (entry.dcDb - entry.hfDb) * p.near_field, e.rate));
            shelf = true;
        } else if (!v.near[ear].identityCoefficients()) {
            v.near[ear].setCoefficients(Biquad{});
            shelf = true;
        }
        if (shelf || !v.near[ear].silentState()) {
            v.near[ear].process(earOut, B);
            e.opReal += 5ull * B;
        }
        applyRamp(earOut, B, firstBlock ? gainTo : v.nearGain[ear], gainTo);
        v.nearGain[ear] = gainTo;
        delayEar(e, v, ear, earOut, v.delayTarget[ear], ear == 0 ? bus.left : bus.right);
    }
    v.delayWrite = (v.delayWrite + B) & (kDelayFrames - 1);
}

// A source on speakers: amplitude panning between the nearest ones, straight into the speaker mix with the gain of its bus chain, as the
// ambisonic mix takes it. Distance and air are already in the signal; head-related cues have no place on speakers.
void renderSpeakers(ms_engine& e, RenderSource& s, const float* mono, float mixFrom, float mixTo, Vec3 direction, std::uint32_t busIndex) {
    const std::uint32_t B = e.B, N = e.speakers;
    float target[kMaxSpeakers];
    e.layout.pan(direction, target);
    if (!s.haveSpeakerGain) {
        std::memcpy(s.speakerGain, target, N * sizeof(float));
        s.haveSpeakerGain = true;
    }
    const float product = e.busProduct[busIndex], productFrom = e.busProductFrom[busIndex];
    const float invB = 1.0f / static_cast<float>(B);
    for (std::uint32_t k = 0; k < N; ++k) {
        const float from = s.speakerGain[k] * mixFrom * productFrom, to = target[k] * mixTo * product;
        s.speakerGain[k] = target[k];
        if (from == 0.0f && to == 0.0f) continue;
        const float step = (to - from) * invB;
        float* out = e.speakerMix.data() + static_cast<std::size_t>(k) * B;
#pragma clang loop vectorize(enable)
        for (std::uint32_t t = 0; t < B; ++t) out[t] += mono[t] * (from + step * static_cast<float>(t + 1));
    }
    e.opReal += static_cast<std::uint64_t>(N) * B;
}

void renderSource(ms_engine& e, std::uint32_t index) {
    RenderSource& s = e.sources[index];
    SourcePublished& out = e.published[index];
    if (s.state == RenderSource::State::Idle) return;
    const std::uint32_t B = e.B;
    const SourceSnapshot& snap = e.frameSources[index];
    const ms_source_params& p = snap.params;

    if (s.scheduled) {
        if (s.startFrame >= e.frame + B) return;
        s.scheduled = false;
        const double offset = s.startFrame > e.frame ? static_cast<double>(s.startFrame - e.frame) : 0.0;
        s.position = s.seekTarget - offset;
        s.seekTarget = 0.0;
    }

    const Vec3 world = Vec3(p.position) - e.listener->position;
    const Vec3 rel{dot(world, e.right), dot(world, e.up), dot(world, e.back)};
    const Vec3 head = p.mode == MS_SPATIAL_PAN ? Vec3{0, 0, -1} : rotate(e.headFrom, rel);
    const float rate = e.rate;
    const float alpha = e.smoothing;
    const bool starting = !s.primed;
    if (starting) {
        s.smoothed = head;
        s.primed = true;
        if (p.propagation_delay) s.position -= static_cast<double>(length(world)) / e.desc.speed_of_sound * rate;
    } else {
        s.smoothed = s.smoothed + (head - s.smoothed) * alpha;
    }
    const float distance = length(s.smoothed);
    const Vec3 direction = normalized(s.smoothed);
    s.smoothedDirection = direction;

    // Pitch: Doppler from the velocities along the line between listener and source.
    float rateTo = std::clamp(p.pitch, kMinPitch, 4.0f);
    if (p.doppler > 0.0f && p.mode != MS_SPATIAL_PAN) {
        const Vec3 toSource = normalized(world);
        const float c = e.desc.speed_of_sound, limit = 0.9f * c;
        const float vs = std::clamp(dot(p.velocity, toSource) * p.doppler, -limit, limit);
        const float vl = std::clamp(dot(e.listener->velocity, toSource) * p.doppler, -limit, limit);
        rateTo *= (c + vl) / (c + vs);
    }
    rateTo = std::clamp(rateTo, kMinRate, 4.0f);
    const float rateFrom = starting ? rateTo : s.rate;
    s.rate = rateTo;

    float* mono = e.scratch.data() + 8ull * B + 512 + 2ull * B;
    if (s.state == RenderSource::State::Tail) std::memset(mono, 0, B * sizeof(float));
    else if (s.callback) readCallback(e, s, mono, rateFrom, rateTo);
    else readBuffer(e, s, mono, rateFrom, rateTo);

    // A callback may hand us anything; one bad sample must not poison the filters.
    if (!finiteBlock(mono, B)) {
        std::memset(mono, 0, B * sizeof(float));
        clearBiquads(s);
        e.stats.faults.fetch_add(1, std::memory_order_relaxed);
    }

    // Sends hear the source itself, as the reflections do. The filters, EQ and air are the path to the listener, so a source
    // behind a door still excites the room it stands in. They do share the source's gain, distance and directivity.
    bool sending = false;
    for (std::uint32_t i = 0; i < MS_MAX_SENDS; ++i) sending = sending || (p.send_bus[i] && p.send_gain[i] != 0.0f);
    float* send = mono + 2ull * B;
    if (sending) std::memcpy(send, mono, B * sizeof(float));

    float* dry = mono + B;  // The signal the reflections hear: before occlusion and distance.
    std::memcpy(dry, mono, B * sizeof(float));
    const float dryTo = s.state == RenderSource::State::Stopping ? 0.0f : p.gain;
    applyRamp(dry, B, starting ? dryTo : s.dryGain, dryTo);
    s.dryGain = dryTo;

    updateStage(s, RenderSource::HighPass, p.highpass_hz, 0.5f, [&] { return p.highpass_hz > 0.0f ? designHighpass(p.highpass_hz, 0.7071f, rate) : Biquad{}; });
    updateStage(s, RenderSource::LowPass, p.lowpass_hz, 0.5f, [&] { return p.lowpass_hz > 0.0f ? designLowpass(p.lowpass_hz, 0.7071f, rate) : Biquad{}; });
    updateStage(s, RenderSource::EqLow, p.eq_low_db, 0.05f, [&] { return p.eq_low_db != 0.0f ? designLowShelf(400.0f, p.eq_low_db, rate) : Biquad{}; });
    updateStage(s, RenderSource::EqMid, p.eq_mid_db, 0.05f, [&] { return p.eq_mid_db != 0.0f ? designPeak(1000.0f, 0.7f, p.eq_mid_db, rate) : Biquad{}; });
    updateStage(s, RenderSource::EqHigh, p.eq_high_db, 0.05f, [&] { return p.eq_high_db != 0.0f ? designHighShelf(2500.0f, p.eq_high_db, rate) : Biquad{}; });

    float directivity = 1.0f, rearAmount = 0.0f;
    if ((p.directivity_alpha > 0.0f || p.directivity_hf_db != 0.0f) && p.mode != MS_SPATIAL_PAN) {
        const float cosine = dot(normalized(Vec3(p.forward)), normalized(Vec3{} - world));
        const float pattern = std::fabs((1.0f - p.directivity_alpha) + p.directivity_alpha * cosine);
        directivity = std::pow(pattern, std::max(p.directivity_sharpness, 0.0f));
        rearAmount = 0.5f * (1.0f - cosine);
    }
    if (p.air_absorption > 0.0f && p.mode != MS_SPATIAL_PAN) {
        const float amount = p.air_absorption * distance;
        const float shelf1 = -e.airTable[0] * amount, shelf2 = -e.airTable[1] * amount;
        updateStage(s, RenderSource::Air1, shelf1, 0.05f, [&] { return shelf1 < -0.02f ? designHighShelf(kAirShelfHz[0], shelf1, rate) : Biquad{}; });
        updateStage(s, RenderSource::Air2, shelf2, 0.05f, [&] { return shelf2 < -0.02f ? designHighShelf(kAirShelfHz[1], shelf2, rate) : Biquad{}; });
    } else {
        updateStage(s, RenderSource::Air1, 0.0f, 0.05f, [&] { return Biquad{}; });
        updateStage(s, RenderSource::Air2, 0.0f, 0.05f, [&] { return Biquad{}; });
    }
    const float rearDb = rearAmount * p.directivity_hf_db;
    updateStage(s, RenderSource::Rear, rearDb, 0.05f, [&] { return rearDb != 0.0f ? designHighShelf(2000.0f, rearDb, rate) : Biquad{}; });
    processStages(e, s, mono);

    const float distanceFactor = p.mode == MS_SPATIAL_PAN ? 1.0f : distanceGain(p, distance);
    float gainTo = std::max(p.gain * distanceFactor * directivity, 0.0f);
    if (s.state == RenderSource::State::Stopping) gainTo = 0.0f;
    const float gainFrom = starting ? gainTo : s.gain;
    applyRamp(mono, B, gainFrom, gainTo);
    if (sending) applyRamp(send, B, gainFrom, gainTo);
    s.gain = gainTo;
    s.score = gainTo;
    e.opReal += B;

    // Everything is checked before it can reach a shared bus, so only this source pays for a bad parameter or filter.
    if (!finiteBlock(mono, B) || !finiteBlock(dry, B) || (sending && !finiteBlock(send, B))) {
        std::memset(mono, 0, 3ull * B * sizeof(float));
        clearBiquads(s);
        s.gain = s.dryGain = 0.0f;
        e.stats.faults.fetch_add(1, std::memory_order_relaxed);
    }

    float level = 0.0f;
    vDSP_maxmgv(mono, 1, &level, B);
    s.peak = level;

    for (std::uint32_t i = 0; i < MS_MAX_SENDS; ++i) {
        const float g = p.send_gain[i];
        if (!p.send_bus[i] || g == 0.0f) continue;
        RenderBus& bus = e.buses[std::min<std::uint32_t>((p.send_bus[i] & 0xffffu) - 1, e.desc.max_buses - 1)];
        vDSP_vsma(send, 1, &g, bus.left, 1, bus.left, 1, B);
        vDSP_vsma(send, 1, &g, bus.right, 1, bus.right, 1, B);
        e.opReal += 2ull * B;
    }

    const std::uint32_t busIndex = std::min<std::uint32_t>(p.bus ? (p.bus & 0xffffu) - 1 : 0, e.desc.max_buses - 1);
    if (p.mode == MS_SPATIAL_PAN) {
        RenderBus& bus = e.buses[busIndex];
        const float angle = (std::clamp(p.pan, -1.0f, 1.0f) + 1.0f) * 0.25f * kPi;
        const float toL = std::cos(angle), toR = std::sin(angle);
        const float fromL = starting ? toL : s.pan[0], fromR = starting ? toR : s.pan[1];
        const float stepL = (toL - fromL) / static_cast<float>(B), stepR = (toR - fromR) / static_cast<float>(B);
        for (std::uint32_t t = 0; t < B; ++t) {
            bus.left[t] += mono[t] * (fromL + stepL * static_cast<float>(t + 1));
            bus.right[t] += mono[t] * (fromR + stepR * static_cast<float>(t + 1));
        }
        s.pan[0] = toL;
        s.pan[1] = toR;
        e.opReal += 2ull * B;
    } else {
        // Direct or shared mix. A voice gained or lost crosses over in one block.
        bool direct = p.mode == MS_SPATIAL_HRTF && s.wantDirect && (e.set || e.speakers);
        if (direct && !s.voice) acquireVoice(e, s);
        direct = direct && s.voice;
        if (direct && starting) s.directMix = 1.0f;
        const float mixFrom = s.directMix, mixTo = direct ? 1.0f : 0.0f;
        if (s.voice) {
            // A demoted voice stops taking input but keeps its state, so its filter tail and every reflection already in the
            // delay lines ring out. The voice is free again after that, and a promotion in the meantime simply resumes it.
            if (mixTo > 0.0f) s.voiceDrain = 0;
            else if (mixFrom > 0.0f) s.voiceDrain = voiceDrainBlocks(e, *s.voice);
            if (e.speakers) renderSpeakers(e, s, mono, mixFrom, mixTo, direction, busIndex);
            else renderDirect(e, s, snap, mono, mixFrom, mixTo, distance, direction);
            applyRamp(dry, B, mixFrom, mixTo);
            renderReflections(e, *s.voice, snap, dry, e.busProduct[busIndex]);
            if (mixTo == 0.0f && (s.voiceDrain == 0 || --s.voiceDrain == 0)) releaseVoice(e, s);
        }
        s.directMix = mixTo;
        const float bedTo = 1.0f - mixTo;
        const float bedFrom = starting ? bedTo : (1.0f - mixFrom) * s.encodeWeight;
        if (e.C > 0 && (bedFrom > 0.0f || bedTo > 0.0f)) {
            float encodeTarget[kMaxChannels];
            const std::uint32_t order = e.desc.ambisonic_order;
            encodeDirection(headToAmbisonic(direction), order, encodeTarget);
            const float product = e.busProduct[busIndex], productFrom = e.busProductFrom[busIndex];
            const float invB = 1.0f / static_cast<float>(B);
            if (s.encodeWeight == 0.0f) std::memcpy(s.encode, encodeTarget, e.C * sizeof(float));  // New to the mix: no pan sweep.
            for (std::uint32_t c = 0; c < e.C; ++c) {
                const float from = s.encode[c] * bedFrom * productFrom, to = encodeTarget[c] * bedTo * product, step = (to - from) * invB;
                float* mix = e.shMix.data() + static_cast<std::size_t>(c) * B;
#pragma clang loop vectorize(enable)
                for (std::uint32_t t = 0; t < B; ++t) mix[t] += mono[t] * (from + step * static_cast<float>(t + 1));
                s.encode[c] = encodeTarget[c];
            }
            e.opReal += static_cast<std::uint64_t>(e.C) * 2 * B;
            s.encodeWeight = bedTo > 0.0f ? 1.0f : 0.0f;
            e.shActiveBlocks = mixHoldBlocks(e);
        }
    }

    // Finished content rings through the filters, then the source goes idle.
    bool ending = s.contentEnded || s.state == RenderSource::State::Stopping;
    if (ending && s.state != RenderSource::State::Tail) {
        s.state = RenderSource::State::Tail;
        bool anyTap = snap.reflectionCount > 0;
        if (s.voice)
            for (const HrtfVoice::Tap& t : s.voice->taps) anyTap = anyTap || t.active;
        s.tailBlocks = e.maxPartitions + (kDelayFrames + B - 1) / B + (anyTap ? kReflectionFrames / B : 0);
        s.contentEnded = true;
    } else if (s.state == RenderSource::State::Tail) {
        if (s.tailBlocks == 0 || --s.tailBlocks == 0) {
            s.state = RenderSource::State::Idle;
            releaseVoice(e, s);
        }
    }
    out.playing.store(s.state != RenderSource::State::Idle, std::memory_order_relaxed);
    out.direct.store(s.voice != nullptr && s.directMix > 0.0f, std::memory_order_relaxed);
    out.position.store(s.position, std::memory_order_relaxed);
    out.level.store(level, std::memory_order_relaxed);
}

void rotationToHead(ms_engine& e, const BedSnapshot& snap, float m[9]) {
    // A bed-frame ambisonic vector -> head-frame ambisonic vector. Conversions: ambisonic (x, y, z) = (-hz, -hx, hy).
    const Quat orientation = normalized(snap.orientation);
    const auto column = [&](Vec3 ambisonicAxis) {
        const Vec3 bedHead{-ambisonicAxis.y, ambisonicAxis.z, -ambisonicAxis.x};
        Vec3 world = rotate(orientation, bedHead);
        Vec3 head = world;
        if (!snap.headLocked) {
            const Vec3 rel{dot(world, e.right), dot(world, e.up), dot(world, e.back)};
            head = rotate(e.headFrom, rel);
        }
        return headToAmbisonic(head);
    };
    const Vec3 cx = column({1, 0, 0}), cy = column({0, 1, 0}), cz = column({0, 0, 1});
    m[0] = cx.x; m[1] = cy.x; m[2] = cz.x;
    m[3] = cx.y; m[4] = cy.y; m[5] = cz.y;
    m[6] = cx.z; m[7] = cy.z; m[8] = cz.z;
}

void renderBed(ms_engine& e, std::uint32_t index) {
    RenderBed& b = e.beds[index];
    if (b.state == RenderBed::State::Idle) return;
    const std::uint32_t B = e.B;
    if (b.scheduled) {
        if (b.startFrame >= e.frame + B) return;
        b.scheduled = false;
        b.position = b.startFrame > e.frame ? -static_cast<std::int64_t>(b.startFrame - e.frame) : 0;
    }
    const BedSnapshot& snap = e.frameBeds[index];
    const std::uint32_t order = std::min(b.order, e.desc.ambisonic_order);
    const std::uint32_t channels = (order + 1) * (order + 1);
    float* in = e.shInput.data();  // [channels][B]
    bool ended = false;
    if (b.callback) {
        std::uint32_t lead = 0;
        if (b.position < 0) {
            lead = std::min<std::uint32_t>(B, static_cast<std::uint32_t>(-b.position));
            b.position += lead;
        }
        float* interleaved = e.shInput.data() + static_cast<std::size_t>(kMaxChannels) * B;
        const std::uint32_t sourceChannels = (b.order + 1) * (b.order + 1);
        std::memset(interleaved, 0, static_cast<std::size_t>(B) * sourceChannels * sizeof(float));
        if (B > lead) {
            const std::uint32_t got = b.callback(b.user, interleaved + static_cast<std::size_t>(lead) * sourceChannels, sourceChannels, B - lead);
            if (got < B - lead) ended = true;
        }
        for (std::uint32_t c = 0; c < channels; ++c) {
            float* target = in + static_cast<std::size_t>(c) * B;
            for (std::uint32_t t = 0; t < B; ++t) target[t] = interleaved[static_cast<std::size_t>(t) * sourceChannels + c];
        }
    } else if (b.buffer) {
        const BufferObject& buffer = *b.buffer;
        const std::uint32_t stride = buffer.channels;
        const auto frames = static_cast<std::int64_t>(buffer.frames);
        for (std::uint32_t c = 0; c < channels; ++c) {
            float* target = in + static_cast<std::size_t>(c) * B;
            const float* data = buffer.data.data() + c;
            std::int64_t position = b.position;
            for (std::uint32_t t = 0; t < B; ++t) {
                if (b.loop && position >= frames) position -= frames;
                target[t] = position >= 0 && position < frames && c < stride ? data[position * stride] : 0.0f;
                ++position;
            }
        }
        b.position += B;
        if (b.loop) { if (b.position >= frames) b.position -= frames; }
        else if (b.position >= frames) ended = true;
    } else {
        ended = true;
    }
    // Content is arbitrary. A bad block must not reach the decoder's history, which would mute the mix for its whole length.
    if (!finiteBlock(in, channels * B)) {
        std::memset(in, 0, static_cast<std::size_t>(channels) * B * sizeof(float));
        e.stats.faults.fetch_add(1, std::memory_order_relaxed);
    }

    float next[kRotationCells], matrix[9];
    rotationToHead(e, snap, matrix);
    e.rotation.compute(matrix, order, next);
    if (!b.haveRotation) {
        std::memcpy(b.rotation, next, sizeof(float) * kRotationCells);
        b.haveRotation = true;
    }
    float gainTo = b.state == RenderBed::State::Stopping ? 0.0f : snap.gain;
    const float gainFrom = b.haveGain ? b.gain : gainTo;
    b.haveGain = true;
    const float invB = 1.0f / static_cast<float>(B);
    const float product = e.busProduct[snap.bus], productFrom = e.busProductFrom[snap.bus];
    std::size_t cell = 0;
    for (std::uint32_t l = 0; l <= order; ++l) {
        const std::uint32_t size = 2 * l + 1;
        for (std::uint32_t i = 0; i < size; ++i) {
            float* mix = e.shMix.data() + static_cast<std::size_t>(l * l + i) * B;
            for (std::uint32_t j = 0; j < size; ++j, ++cell) {
                const float from = b.rotation[cell] * gainFrom * productFrom, to = next[cell] * gainTo * product;
                if (std::fabs(from) < 1e-7f && std::fabs(to) < 1e-7f) continue;
                const float step = (to - from) * invB;
                const float* source = in + static_cast<std::size_t>(l * l + j) * B;
#pragma clang loop vectorize(enable)
                for (std::uint32_t t = 0; t < B; ++t) mix[t] += source[t] * (from + step * static_cast<float>(t + 1));
                e.opReal += B;
            }
        }
    }
    std::memcpy(b.rotation, next, cell * sizeof(float));
    b.gain = gainTo;
    e.shActiveBlocks = mixHoldBlocks(e);
    if (ended || b.state == RenderBed::State::Stopping) b.state = RenderBed::State::Idle;
}

void decodeMix(ms_engine& e) {
    if (e.speakers) {
        // AllRAD: one matrix row per speaker, no filters and no latency.
        if (e.shActiveBlocks == 0 || e.C == 0) return;
        --e.shActiveBlocks;
        const std::uint32_t B = e.B, C = std::min(e.C, e.layout.decodeChannels());
        for (std::uint32_t k = 0; k < e.speakers; ++k) {
            const float* row = e.layout.decodeRow(k);
            float* out = e.speakerMix.data() + static_cast<std::size_t>(k) * B;
            for (std::uint32_t c = 0; c < C; ++c)
                if (row[c] != 0.0f) vDSP_vsma(e.shMix.data() + static_cast<std::size_t>(c) * B, 1, &row[c], out, 1, out, 1, B);
        }
        e.opReal += static_cast<std::uint64_t>(e.speakers) * C * B;
        return;
    }
    if (!e.decoder || e.shActiveBlocks == 0) return;
    --e.shActiveBlocks;
    const std::uint32_t B = e.B, C = e.C;
    const AmbisonicDecoder& decoder = *e.decoder;
    e.shFdlHead = (e.shFdlHead + 1) % decoder.partitions;
    for (std::uint32_t c = 0; c < C; ++c) {
        float* window = e.shInput.data() + static_cast<std::size_t>(kMaxChannels) * B * 2 + static_cast<std::size_t>(c) * 2 * B;
        std::memcpy(window, window + B, B * sizeof(float));
        std::memcpy(window + B, e.shMix.data() + static_cast<std::size_t>(c) * B, B * sizeof(float));
        const std::size_t at = (static_cast<std::size_t>(c) * decoder.partitions + e.shFdlHead) * B;
        e.fft.forward(window, e.shFdlRe.data() + at, e.shFdlIm.data() + at);
    }
    e.opForward += C;
    float* accRe = e.scratch.data() + 2ull * B;
    float* accIm = accRe + B;
    float* time = e.scratch.data() + 4ull * B;
    float* older = e.scratch.data() + 6ull * B;
    const auto accumulate = [&](const AmbisonicDecoder& d, std::uint32_t ear) {
        std::memset(accRe, 0, 2ull * B * sizeof(float));
        for (std::uint32_t c = 0; c < std::min(C, d.channels); ++c)
            for (std::uint32_t p = 0; p < d.partitions; ++p) {
                const std::uint32_t age = (e.shFdlHead + decoder.partitions - std::min(p, decoder.partitions - 1)) % decoder.partitions;
                const std::size_t at = (static_cast<std::size_t>(c) * decoder.partitions + age) * B;
                cmacPacked(accRe, accIm, e.shFdlRe.data() + at, e.shFdlIm.data() + at, d.spectrumRe(ear, c, p), d.spectrumIm(ear, c, p), B);
            }
        e.opSpectrum += static_cast<std::uint64_t>(C) * d.partitions * B;
    };
    float* out[2] = {e.buses[0].left, e.buses[0].right};
    for (std::uint32_t ear = 0; ear < 2; ++ear) {
        accumulate(decoder, ear);
        e.fft.inverse(accRe, accIm, time);
        ++e.opInverse;
        if (e.previousDecoder) {
            accumulate(*e.previousDecoder, ear);
            e.fft.inverse(accRe, accIm, older);
            ++e.opInverse;
            const float* ramp = e.crossfadeRamp.data();
            for (std::uint32_t t = 0; t < B; ++t) out[ear][t] += e.inverseScale * (time[B + t] * ramp[t] + older[B + t] * (1.0f - ramp[t]));
        } else {
            vDSP_vsma(time + B, 1, &e.inverseScale, out[ear], 1, out[ear], 1, B);
        }
    }
}

void rankVoices(ms_engine& e) {
    // Keep the loudest sources in direct convolution. Stable, so the ranking is a pure function of the state.
    std::uint32_t count = 0;
    for (std::uint32_t i = 0; i < e.desc.max_sources; ++i) {
        RenderSource& s = e.sources[i];
        s.wantDirect = false;
        if (s.state == RenderSource::State::Idle || e.frameSources[i].params.mode != MS_SPATIAL_HRTF) continue;
        e.ranking[count++] = i;
        if (!s.primed) {  // Nothing rendered yet to score it by, so what it will be: a new sound must not take a voice for free.
            const ms_source_params& p = e.frameSources[i].params;
            s.score = std::max(p.gain * distanceGain(p, length(Vec3(p.position) - e.listener->position)), 0.0f);
        }
    }
    const std::uint32_t keep = std::min(count, e.desc.max_hrtf_voices);
    const RenderSource* sources = e.sources.data();
    const auto rankOf = [sources](std::uint32_t i) {
        const RenderSource& s = sources[i];
        return s.voice && s.directMix > 0.0f ? s.score * kVoiceHoldMargin : s.score;
    };
    std::partial_sort(e.ranking.data(), e.ranking.data() + keep, e.ranking.data() + count, [rankOf](std::uint32_t a, std::uint32_t b) {
        const float scoreA = rankOf(a), scoreB = rankOf(b);
        return scoreA != scoreB ? scoreA > scoreB : a < b;
    });
    for (std::uint32_t i = 0; i < keep; ++i) e.sources[e.ranking[i]].wantDirect = true;
}

void mixBuses(ms_engine& e) {
    const std::uint32_t B = e.B;
    for (std::uint32_t i = e.desc.max_buses; i-- > 0;) {
        RenderBus& bus = e.buses[i];
        if (!bus.used) continue;
        const BusSnapshot& snap = e.frameBuses[i];
        if (snap.highpass != bus.highpassHz) {
            bus.highpassHz = snap.highpass;
            const Biquad q = snap.highpass > 0.0f ? designHighpass(snap.highpass, 0.7071f, e.rate) : Biquad{};
            bus.highpass[0].setCoefficients(q);
            bus.highpass[1].setCoefficients(q);
        }
        if (snap.lowpass != bus.lowpassHz) {
            bus.lowpassHz = snap.lowpass;
            const Biquad q = snap.lowpass > 0.0f ? designLowpass(snap.lowpass, 0.7071f, e.rate) : Biquad{};
            bus.lowpass[0].setCoefficients(q);
            bus.lowpass[1].setCoefficients(q);
        }
        float* channels[2] = {bus.left, bus.right};
        for (int c = 0; c < 2; ++c) {
            if (!bus.highpass[c].identityCoefficients()) bus.highpass[c].process(channels[c], B);
            if (!bus.lowpass[c].identityCoefficients()) bus.lowpass[c].process(channels[c], B);
            applyRamp(channels[c], B, bus.primed ? bus.gain : snap.gain, snap.gain);
        }
        bus.gain = snap.gain;
        bus.primed = true;
        if (bus.processor) bus.processor(bus.user, bus.left, bus.right, B, e.frame);
        // The last gate before a shared sum: a bad child bus loses its own block, and the master keeps everything else.
        if (i != 0 && (!finiteBlock(bus.left, B) || !finiteBlock(bus.right, B))) {
            std::memset(bus.left, 0, B * sizeof(float));
            std::memset(bus.right, 0, B * sizeof(float));
            for (int c = 0; c < 2; ++c) {
                bus.highpass[c].clear();
                bus.lowpass[c].clear();
            }
            e.stats.faults.fetch_add(1, std::memory_order_relaxed);
        }
        if (i != 0) {
            RenderBus& parent = e.buses[bus.parent];
            vDSP_vadd(parent.left, 1, bus.left, 1, parent.left, 1, B);
            vDSP_vadd(parent.right, 1, bus.right, 1, parent.right, 1, B);
        }
    }
}

void renderBlock(ms_engine& e, float* const* outs, std::uint32_t outCount) {
    const std::uint32_t B = e.B;
    processCommands(e);  // Commands before the snapshot: any command read has its snapshot already published.
    if (e.blockCallback) e.blockCallback(e.blockUser, B, e.frame);
    e.triple.acquire();
    const std::uint32_t r = e.triple.readIndex();
    e.listener = &e.snapListener[r];
    e.frameSources = e.snapSources[r].data();
    e.frameBeds = e.snapBeds[r].data();
    e.frameBuses = e.snapBuses[r].data();

    const ListenerSnapshot& l = *e.listener;
    const Vec3 forward = normalized(l.forward);
    e.right = normalized(cross(forward, l.up), {1, 0, 0});
    e.up = cross(e.right, forward);
    e.back = forward * -1.0f;
    e.headFrom = conjugate(unpackQuat(e.headPacked.load(std::memory_order_acquire)));
    if (e.speakers) {
        e.headFrom = Quat{};  // Speakers stand still in the room whichever way the head turns.
        std::memset(e.speakerMix.data(), 0, static_cast<std::size_t>(e.speakers) * B * sizeof(float));
    }

    for (std::uint32_t i = 0; i < e.desc.max_buses; ++i) {
        if (!e.buses[i].used) continue;
        std::memset(e.buses[i].left, 0, B * sizeof(float));
        std::memset(e.buses[i].right, 0, B * sizeof(float));
        const float product = i == 0 ? 1.0f : e.frameBuses[i].gain * e.busProduct[e.buses[i].parent];
        e.busProductFrom[i] = e.buses[i].primed ? e.busProduct[i] : product;
        e.busProduct[i] = product;
    }
    if (e.C) std::memset(e.shMix.data(), 0, static_cast<std::size_t>(e.C) * B * sizeof(float));

    if ((e.frame / B) % kRankInterval == 0) rankVoices(e);
    else {
        // A source that starts between rankings takes a free voice straight away.
        for (std::uint32_t i = 0; i < e.desc.max_sources; ++i)
            if (e.sources[i].state != RenderSource::State::Idle && !e.sources[i].primed && e.freeVoiceCount > 0) e.sources[i].wantDirect = true;
    }
    for (std::uint32_t i = 0; i < e.desc.max_beds; ++i) renderBed(e, i);
    std::uint32_t playing = 0, direct = 0, shared = 0;
    for (std::uint32_t i = 0; i < e.desc.max_sources; ++i) {
        renderSource(e, i);
        if (e.sources[i].state != RenderSource::State::Idle) {
            ++playing;
            if (e.sources[i].voice && e.sources[i].directMix > 0.0f) ++direct;
            else ++shared;
        }
    }
    decodeMix(e);
    e.previousDecoder = nullptr;
    if (e.retiringSet) {
        if (e.retired.push({Retired::Kind::Hrtf, e.retiringSet, nullptr})) e.retiringSet = nullptr;
    }
    const float masterTo = e.frameBuses[0].gain, masterFrom = e.buses[0].primed ? e.buses[0].gain : masterTo;
    mixBuses(e);

    RenderBus& master = e.buses[0];
    // The finished channels: the master's two for headphones, or the speakers, which take the master's gain here and its stereo content
    // as a phantom pair.
    float* finals[kMaxSpeakers];
    std::uint32_t F = 2;
    if (e.speakers) {
        F = e.speakers;
        for (std::uint32_t k = 0; k < F; ++k) {
            float* out = finals[k] = e.speakerMix.data() + static_cast<std::size_t>(k) * B;
            applyRamp(out, B, masterFrom, masterTo);
            const float gl = e.layout.stereoIn(0)[k], gr = e.layout.stereoIn(1)[k];
            if (gl != 0.0f) vDSP_vsma(master.left, 1, &gl, out, 1, out, 1, B);
            if (gr != 0.0f) vDSP_vsma(master.right, 1, &gr, out, 1, out, 1, B);
        }
        e.opReal += 3ull * F * B;
    } else {
        finals[0] = master.left;
        finals[1] = master.right;
    }
    if (e.desc.limiter) {
        float peak = 0.0f;
        for (std::uint32_t c = 0; c < F; ++c) {
            float channelPeak = 0.0f;
            vDSP_maxmgv(finals[c], 1, &channelPeak, B);
            peak = std::max(peak, channelPeak);
        }
        const float target = peak > 1.0f ? 1.0f / peak : 1.0f;
        const float next = target < e.limiterGain ? target : std::min(1.0f, e.limiterGain + 0.002f);
        const float low = -1.0f, high = 1.0f;
        for (std::uint32_t c = 0; c < F; ++c) {
            applyRamp(finals[c], B, e.limiterGain, next);
            vDSP_vclip(finals[c], 1, &low, &high, finals[c], 1, B);
        }
        e.limiterGain = next;
    }
    float sums[2] = {0, 0}, peaks[2] = {0, 0};
    bool finite = true;
    std::uint64_t overs = 0;
    for (std::uint32_t c = 0; c < F; ++c) {
        float sum = 0.0f, peak = 0.0f;
        vDSP_svesq(finals[c], 1, &sum, B);
        vDSP_maxmgv(finals[c], 1, &peak, B);
        finite = finite && std::isfinite(sum);
        if (c < 2) {
            sums[c] = sum;
            peaks[c] = peak;
        }
        if (peak > 1.0f)
            for (std::uint32_t t = 0; t < B; ++t) overs += std::fabs(finals[c][t]) > 1.0f;
    }
    if (!finite) {
        for (std::uint32_t c = 0; c < F; ++c) std::memset(finals[c], 0, B * sizeof(float));
        e.stats.faults.fetch_add(1, std::memory_order_relaxed);
        sums[0] = sums[1] = peaks[0] = peaks[1] = 0.0f;
        overs = 0;
    }
    if (F == 1) {
        sums[1] = sums[0];
        peaks[1] = peaks[0];
    }
    if (outCount == F) {
        for (std::uint32_t c = 0; c < F; ++c) std::memcpy(outs[c], finals[c], B * sizeof(float));
    } else {
        // A speaker render asked for as stereo: fold the speakers down.
        for (std::uint32_t side = 0; side < 2; ++side) {
            std::memset(outs[side], 0, B * sizeof(float));
            for (std::uint32_t k = 0; k < F; ++k) {
                const float g = e.layout.stereoOut(side)[k];
                if (g != 0.0f) vDSP_vsma(finals[k], 1, &g, outs[side], 1, outs[side], 1, B);
            }
        }
    }
    e.frame += B;

    AtomicStats& st = e.stats;
    const auto relaxed = std::memory_order_relaxed;
    st.frames.store(e.frame, relaxed);
    st.blocks.fetch_add(1, relaxed);
    st.overs.fetch_add(overs, relaxed);
    st.playing.store(playing, relaxed);
    st.direct.store(direct, relaxed);
    st.ambisonic.store(shared, relaxed);
    for (int c = 0; c < 2; ++c) {
        st.peak[c].store(peaks[c], relaxed);
        st.rms[c].store(std::sqrt(sums[c] / static_cast<float>(B)), relaxed);
    }
    st.opBlocks.fetch_add(1, relaxed);
    st.forward.store(e.opForward, relaxed);
    st.inverse.store(e.opInverse, relaxed);
    st.spectrum.store(e.opSpectrum, relaxed);
    st.real.store(e.opReal, relaxed);
    st.updates.store(e.opUpdates, relaxed);
}

}  // namespace

namespace {

// Renders frames into outCount channels: whole blocks straight into the caller's memory, a partial one through the FIFO.
std::uint64_t renderFrames(ms_engine& e, float* const* outs, std::uint32_t outCount, std::uint32_t frames) {
    const std::uint64_t begin = mach_absolute_time();
    const std::uint32_t B = e.B;
    float* fifo[kMaxSpeakers];
    for (std::uint32_t c = 0; c < outCount; ++c) fifo[c] = e.fifo.data() + static_cast<std::size_t>(c) * B;
    if (e.fifoChannels != outCount) e.fifoRead = e.fifoFill;  // Left by the other entry point, in another channel count.
    e.fifoChannels = outCount;
    float* at[kMaxSpeakers];
    std::uint32_t done = 0;
    while (done < frames) {
        if (e.fifoRead < e.fifoFill) {
            const std::uint32_t take = std::min(frames - done, e.fifoFill - e.fifoRead);
            for (std::uint32_t c = 0; c < outCount; ++c) std::memcpy(outs[c] + done, fifo[c] + e.fifoRead, take * sizeof(float));
            e.fifoRead += take;
            done += take;
        } else if (frames - done >= B) {
            for (std::uint32_t c = 0; c < outCount; ++c) at[c] = outs[c] + done;
            renderBlock(e, at, outCount);
            done += B;
        } else {
            renderBlock(e, fifo, outCount);
            e.fifoRead = 0;
            e.fifoFill = B;
        }
    }
    const double seconds = static_cast<double>(mach_absolute_time() - begin) * e.ticksToSeconds;
    e.stats.last.store(seconds, std::memory_order_relaxed);
    if (seconds > e.stats.max.load(std::memory_order_relaxed)) e.stats.max.store(seconds, std::memory_order_relaxed);
    return e.frame;
}

}  // namespace

extern "C" uint64_t ms_engine_render(ms_engine* engine, float* left, float* right, uint32_t frames) {
    float* outs[2] = {left, right};
    return renderFrames(*engine, outs, 2, frames);
}

extern "C" uint64_t ms_engine_render_speakers(ms_engine* engine, float* const* channels, uint32_t frames) {
    return renderFrames(*engine, channels, engine->speakers ? engine->speakers : 2, frames);
}

extern "C" uint32_t ms_engine_output_channels(const ms_engine* engine) { return engine->speakers ? engine->speakers : 2; }

extern "C" uint64_t ms_engine_render_interleaved(ms_engine* engine, float* stereo, uint32_t frames) {
    float left[1024], right[1024];
    uint64_t position = engine->frame;
    for (std::uint32_t done = 0; done < frames;) {
        const std::uint32_t chunk = std::min<std::uint32_t>(frames - done, 1024);
        position = ms_engine_render(engine, left, right, chunk);
        const DSPSplitComplex split{left, right};
        vDSP_ztoc(&split, 1, reinterpret_cast<DSPComplex*>(stereo + 2ull * done), 2, chunk);
        done += chunk;
    }
    return position;
}
