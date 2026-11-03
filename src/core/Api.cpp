// SPDX-License-Identifier: Apache-2.0
// The control half: everything in the public C API except rendering.
#include <AudioToolbox/AudioToolbox.h>
#include <mach/mach_time.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <string>

#include "Engine.hpp"

using namespace mspatial::core;
using mspatial::hrtf::HrtfSet;

namespace {

constexpr std::uint32_t kIndexMask = 0xffffu;

struct Guard {
    os_unfair_lock* lock;
    explicit Guard(ms_engine* e) : lock(&e->lock) { os_unfair_lock_lock(lock); }
    ~Guard() { os_unfair_lock_unlock(lock); }
};

inline HrtfSet* asSet(ms_hrtf* h) { return reinterpret_cast<HrtfSet*>(h); }
inline const HrtfSet* asSet(const ms_hrtf* h) { return reinterpret_cast<const HrtfSet*>(h); }
inline ms_hrtf* asHandle(HrtfSet* s) { return reinterpret_cast<ms_hrtf*>(s); }

std::uint32_t allocateSlot(std::vector<ms_engine::Slot>& slots, std::uint32_t first = 0) {
    for (std::uint32_t i = first; i < slots.size(); ++i)
        if (!slots[i].alive) {
            slots[i].alive = true;
            return (static_cast<std::uint32_t>(slots[i].generation) << 16) | (i + 1);
        }
    return 0;
}

bool validSlot(const std::vector<ms_engine::Slot>& slots, std::uint32_t id, std::uint32_t& index) {
    if (id == 0) return false;
    index = (id & kIndexMask) - 1;
    return index < slots.size() && slots[index].alive && slots[index].generation == (id >> 16);
}

void freeSlot(std::vector<ms_engine::Slot>& slots, std::uint32_t index) {
    slots[index].alive = false;
    slots[index].generation = static_cast<std::uint16_t>(slots[index].generation + 1);
    if (slots[index].generation == 0) slots[index].generation = 1;
}

bool finite(ms_vec3 v) { return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z); }

// Everything the render side computes with. A source parameter that is not finite never gets past the API.
bool validParams(const ms_source_params& p) {
    const float values[] = {p.gain, p.pitch, p.pan, p.ref_distance, p.max_distance, p.rolloff, p.doppler, p.air_absorption, p.near_field,
                            p.directivity_alpha, p.directivity_sharpness, p.directivity_hf_db, p.highpass_hz, p.lowpass_hz,
                            p.eq_low_db, p.eq_mid_db, p.eq_high_db};
    for (const float v : values)
        if (!std::isfinite(v)) return false;
    for (const float g : p.send_gain)
        if (!std::isfinite(g)) return false;
    return finite(p.position) && finite(p.velocity) && finite(p.forward) && p.highpass_hz >= 0.0f && p.lowpass_hz >= 0.0f &&
           p.mode <= MS_SPATIAL_PAN && p.distance_model <= MS_DISTANCE_EXPONENTIAL;
}

void collectRetired(ms_engine& e) {
    Retired r;
    while (e.retired.pop(r)) {
        if (r.kind == Retired::Kind::Buffer) {
            for (std::size_t i = 0; i < e.retiringBuffers.size(); ++i)
                if (e.retiringBuffers[i].get() == r.pointer) {
                    e.retiringBuffers.erase(e.retiringBuffers.begin() + static_cast<std::ptrdiff_t>(i));
                    break;
                }
        } else if (r.kind == Retired::Kind::Crosstalk) {
            for (std::size_t i = 0; i < e.heldCrosstalk.size(); ++i)
                if (e.heldCrosstalk[i].get() == r.pointer) {
                    e.heldCrosstalk.erase(e.heldCrosstalk.begin() + static_cast<std::ptrdiff_t>(i));
                    break;
                }
        } else {
            auto* set = static_cast<HrtfSet*>(r.pointer);
            for (std::size_t i = 0; i < e.heldSets.size(); ++i)
                if (e.heldSets[i].set == set) {
                    ms_hrtf_release(asHandle(set));
                    e.heldSets.erase(e.heldSets.begin() + static_cast<std::ptrdiff_t>(i));
                    break;
                }
        }
    }
}

void queue(ms_engine& e, CommandType type, std::uint32_t index, std::uint32_t a = 0, std::uint32_t b = 0, std::uint64_t frame = 0, double position = 0.0,
           void* pointer = nullptr, void* pointer2 = nullptr) {
    e.pending.push_back({type, index, a, b, frame, position, pointer, pointer2});
}

ms_status selectHrtf(ms_engine& e, HrtfSet* set) {
    std::unique_ptr<HrtfSet> converted;
    HrtfSet* chosen = set;
    if (set->sampleRate != e.rate) {
        converted = mspatial::hrtf::resampleSet(*set, e.rate);
        chosen = converted.get();
    } else {
        ms_hrtf_retain(asHandle(set));
    }
    std::shared_ptr<const AmbisonicDecoder> decoder;
    if (e.C) {
        decoder = chosen->decoder(e.desc.ambisonic_order, e.B);
        if (!decoder) {
            if (converted) return MS_ERR_UNSUPPORTED;
            ms_hrtf_release(asHandle(set));
            return MS_ERR_UNSUPPORTED;
        }
    }
    converted.release();
    e.heldSets.push_back({chosen, decoder});
    e.currentSet = chosen;
    e.currentDecoder = decoder;
    queue(e, CommandType::SetHrtf, 0, 0, 0, 0, 0.0, chosen, const_cast<AmbisonicDecoder*>(decoder.get()));
    return MS_OK;
}

// A least-squares fit of two shelves to the ISO curve, in the linear range of the shelves' dB response.
void fitAirShelves(ms_engine& e) {
    constexpr int kProbes = 10;
    double aa = 0, ab = 0, bb = 0, ay = 0, by = 0;
    const auto shelfDb = [&](float corner, double hz) {
        const Biquad q = designHighShelf(corner, 6.0f, e.rate);
        const double w = 2.0 * M_PI * hz / e.rate;
        const std::complex<double> z1 = std::polar(1.0, -w), z2 = z1 * z1;
        const std::complex<double> h = (static_cast<double>(q.b0) + static_cast<double>(q.b1) * z1 + static_cast<double>(q.b2) * z2) /
                                       (1.0 + static_cast<double>(q.a1) * z1 + static_cast<double>(q.a2) * z2);
        return 20.0 * std::log10(std::abs(h)) / 6.0;
    };
    for (int i = 0; i < kProbes; ++i) {
        const double hz = 1000.0 * std::pow(std::min(16000.0, e.rate * 0.45) / 1000.0, i / (kProbes - 1.0));
        const double target = airAttenuationDbPerMetre(hz, e.desc.air_temperature_c, e.desc.air_humidity_percent);
        const double a = shelfDb(kAirShelfHz[0], hz), b = shelfDb(kAirShelfHz[1], hz);
        aa += a * a; ab += a * b; bb += b * b; ay += a * target; by += b * target;
    }
    const double det = aa * bb - ab * ab;
    e.airTable[0] = static_cast<float>(std::max(0.0, (ay * bb - by * ab) / det));
    e.airTable[1] = static_cast<float>(std::max(0.0, (by * aa - ay * ab) / det));
}

}  // namespace

extern "C" {

const char* ms_version_string(void) { return "0.1.0"; }

const char* ms_status_string(ms_status status) {
    switch (status) {
        case MS_OK: return "ok";
        case MS_ERR_INVALID_ARGUMENT: return "invalid argument";
        case MS_ERR_OUT_OF_MEMORY: return "out of memory";
        case MS_ERR_NOT_FOUND: return "not found";
        case MS_ERR_BAD_FILE: return "bad file";
        case MS_ERR_LIMIT_REACHED: return "limit reached";
        case MS_ERR_DEVICE: return "audio device error";
        case MS_ERR_UNSUPPORTED: return "unsupported";
        case MS_ERR_RATE_MISMATCH: return "device and engine sample rates differ";
    }
    return "unknown";
}

void ms_crosstalk_params_default(ms_crosstalk_params* p) {
    if (!p) return;
    *p = ms_crosstalk_params{};
    p->enabled = 1;
    p->speaker_span_degrees = 15.0f;
    p->head_radius_m = 0.0875f;
    p->max_boost_db = 15.0f;
    p->low_hz = 200.0f;
    p->high_hz = 7000.0f;
}

ms_status ms_engine_set_crosstalk(ms_engine* engine, const ms_crosstalk_params* params) {
    if (!engine) return MS_ERR_INVALID_ARGUMENT;
    ms_engine& e = *engine;
    if (e.speakers) return params && params->enabled ? MS_ERR_UNSUPPORTED : MS_OK;
    // Designed before the lock: it takes a few milliseconds, and nothing it reads changes after creation.
    std::unique_ptr<mspatial::core::CrosstalkFilter> filter;
    if (params && params->enabled) {
        filter = mspatial::core::designCrosstalk(*params, e.rate, e.B, e.desc.speed_of_sound);
        if (!filter || filter->partitions > e.crosstalkPartitions) return MS_ERR_INVALID_ARGUMENT;
    }
    Guard guard(engine);
    if (!filter && !e.currentCrosstalk) return MS_OK;
    e.currentCrosstalk = filter.get();
    if (filter) e.heldCrosstalk.push_back(std::move(filter));
    queue(e, CommandType::SetCrosstalk, 0, 0, 0, 0, 0.0, const_cast<mspatial::core::CrosstalkFilter*>(e.currentCrosstalk));
    return MS_OK;
}

uint32_t ms_engine_crosstalk_latency_frames(const ms_engine* engine) {
    return engine && engine->currentCrosstalk ? engine->currentCrosstalk->latencyFrames : 0;
}

ms_status ms_speaker_layout_preset(ms_speaker_preset preset, ms_speaker_layout* layout) {
    if (!layout) return MS_ERR_INVALID_ARGUMENT;
    *layout = ms_speaker_layout{};
    const auto add = [&](float azimuth, float elevation, bool lfe = false) {
        layout->speakers[layout->count++] = {azimuth, elevation, lfe ? 1u : 0u};
    };
    const auto front = [&] {
        add(30, 0);
        add(-30, 0);
        add(0, 0);
        add(0, 0, true);
    };
    switch (preset) {
    case MS_SPEAKERS_MONO: add(0, 0); break;
    case MS_SPEAKERS_STEREO: add(30, 0); add(-30, 0); break;
    case MS_SPEAKERS_QUAD: add(45, 0); add(-45, 0); add(135, 0); add(-135, 0); break;
    case MS_SPEAKERS_5_1: front(); add(110, 0); add(-110, 0); break;
    case MS_SPEAKERS_7_1: front(); add(90, 0); add(-90, 0); add(135, 0); add(-135, 0); break;
    case MS_SPEAKERS_5_1_4: front(); add(110, 0); add(-110, 0); add(30, 30); add(-30, 30); add(110, 30); add(-110, 30); break;
    case MS_SPEAKERS_7_1_4: front(); add(90, 0); add(-90, 0); add(135, 0); add(-135, 0); add(45, 30); add(-45, 30); add(135, 30); add(-135, 30); break;
    default: return MS_ERR_INVALID_ARGUMENT;
    }
    return MS_OK;
}

ms_status ms_speaker_layout_pan(const ms_speaker_layout* layout, ms_vec3 direction, float* gains) {
    if (!layout || !gains) return MS_ERR_INVALID_ARGUMENT;
    mspatial::core::SpeakerLayout built;
    if (!built.build(*layout, 0)) return MS_ERR_INVALID_ARGUMENT;
    // Listener axes to head axes: the same at the identity listener, with -Z forward.
    built.pan(mspatial::core::Vec3(direction), gains);
    return MS_OK;
}

void ms_engine_desc_default(ms_engine_desc* d) {
    *d = ms_engine_desc{};
    d->api_version = MS_API_VERSION;
    d->sample_rate = 48000.0f;
    d->block_frames = 128;
    d->max_sources = 256;
    d->max_hrtf_voices = 64;
    d->max_buses = 16;
    d->max_beds = 8;
    d->max_buffers = 4096;
    d->ambisonic_order = 3;
    d->speed_of_sound = 343.0f;
    d->smoothing_seconds = 0.005f;
    d->air_temperature_c = 20.0f;
    d->air_humidity_percent = 50.0f;
    d->limiter = 0;
}

ms_status ms_engine_create(const ms_engine_desc* desc, ms_engine** engine) {
    if (!desc || !engine || desc->api_version != MS_API_VERSION) return MS_ERR_INVALID_ARGUMENT;
    const ms_engine_desc& d = *desc;
    if (!(d.sample_rate >= 8000.0f && d.sample_rate <= 192000.0f)) return MS_ERR_INVALID_ARGUMENT;
    if (d.block_frames != 64 && d.block_frames != 128 && d.block_frames != 256 && d.block_frames != 512) return MS_ERR_INVALID_ARGUMENT;
    if (d.max_sources == 0 || d.max_sources > 8192 || d.max_buses == 0 || d.max_buses > 256 || d.max_buffers == 0 || d.max_beds > 256 ||
        d.ambisonic_order > kMaxOrder || !(d.speed_of_sound > 1.0f) || d.smoothing_seconds < 0.0f)
        return MS_ERR_INVALID_ARGUMENT;
    try {
        mach_timebase_info_data_t timebase;
        mach_timebase_info(&timebase);
        auto owned = std::make_unique<ms_engine>();
        ms_engine& e = *owned;
        e.desc = d;
        e.ticksToSeconds = static_cast<double>(timebase.numer) / timebase.denom * 1e-9;
        e.desc.max_hrtf_voices = std::min(d.max_hrtf_voices, d.max_sources);
        e.rate = d.sample_rate;
        e.B = d.block_frames;
        const std::uint32_t B = e.B;
        e.C = (d.ambisonic_order + 1) * (d.ambisonic_order + 1);
        if (d.ambisonic_order == 0) e.C = 0;
        e.maxPartitions = (kMaxFilterTaps + B - 1) / B;
        e.inverseScale = 1.0f / (8.0f * static_cast<float>(B));
        e.fft.init(2 * B);

        e.sourceSlots.resize(d.max_sources);
        e.bedSlots.resize(d.max_beds);
        e.bufferSlots.resize(d.max_buffers);
        e.busSlots.resize(d.max_buses);
        e.busSlots[0].alive = true;
        e.busSlots[0].generation = 0;  // MS_BUS_MASTER is the fixed handle 1.
        e.buffers.resize(d.max_buffers);
        e.workSources.resize(d.max_sources);
        e.workBeds.resize(d.max_beds);
        e.workBuses.resize(d.max_buses);
        for (auto& source : e.workSources) {
            ms_source_params_default(&source.params);
            source.reflectionCount = 0;
        }
        for (auto& bed : e.workBeds) bed = {Quat{}, 1.0f, 0, 0};
        for (auto& bus : e.workBuses) bus = {1.0f, 0.0f, 0.0f};
        for (int i = 0; i < 3; ++i) {
            e.snapSources[i].reset(d.max_sources);
            e.snapBeds[i].reset(std::max(d.max_beds, 1u));
            e.snapBuses[i].reset(d.max_buses);
            std::copy(e.workSources.begin(), e.workSources.end(), e.snapSources[i].data());
            std::copy(e.workBeds.begin(), e.workBeds.end(), e.snapBeds[i].data());
            std::copy(e.workBuses.begin(), e.workBuses.end(), e.snapBuses[i].data());
        }
        e.commands.init(8192);
        e.retired.init(d.max_buffers + 64);
        e.published.reset(d.max_sources);

        e.sources.reset(d.max_sources);
        e.beds.reset(std::max(d.max_beds, 1u));
        e.buses.reset(d.max_buses);
        e.ranking.reset(d.max_sources);
        e.scratch.reset(16ull * B + kMaxFilterTaps + 64);
        if (d.speakers.count > 0) {
            if (!e.layout.build(d.speakers, d.ambisonic_order)) return MS_ERR_INVALID_ARGUMENT;
            e.speakers = e.layout.count();
            e.speakerMix.reset(static_cast<std::size_t>(e.speakers) * B);
        }
        e.fifo.reset(static_cast<std::size_t>(std::max(e.speakers, 2u)) * B);
        e.crosstalkPartitions = (kMaxCrosstalkTaps + B - 1) / B;
        e.crosstalkWindow.reset(4ull * B);
        e.crosstalkFdlRe.reset(2ull * e.crosstalkPartitions * B);
        e.crosstalkFdlIm.reset(2ull * e.crosstalkPartitions * B);
        e.crosstalkOut.reset(6ull * B);
        e.crossfadeRamp.reset(B);
        for (std::uint32_t t = 0; t < B; ++t) {
            const float s = std::sin(0.5f * kPi * (static_cast<float>(t) + 0.5f) / static_cast<float>(B));
            e.crossfadeRamp[t] = s * s;
        }
        e.pullMemory.reset(static_cast<std::size_t>(d.max_sources) * (5ull * B + 32));
        for (std::uint32_t i = 0; i < d.max_sources; ++i) {
            RenderSource& s = e.sources[i];
            s = RenderSource{};
            s.pull = e.pullMemory.data() + static_cast<std::size_t>(i) * (5ull * B + 32);
            s.voice = nullptr;
            for (float& v : s.stageValue) v = NAN;
        }
        e.busMemory.reset(static_cast<std::size_t>(d.max_buses) * 2 * B);
        for (std::uint32_t i = 0; i < d.max_buses; ++i) {
            RenderBus& b = e.buses[i];
            b = RenderBus{};
            b.left = e.busMemory.data() + static_cast<std::size_t>(i) * 2 * B;
            b.right = b.left + B;
            b.gain = 1.0f;
        }
        e.buses[0].used = true;
        e.busProduct[0] = e.busProductFrom[0] = 1.0f;

        const std::size_t P = e.maxPartitions;
        const std::size_t perVoice = 2ull * B + 2 * P * B + 8 * P * B + 4ull * kDelayFrames + 3ull * kReflectionFrames + 256;
        e.voices.reset(e.desc.max_hrtf_voices);
        e.voiceMemory.reset(perVoice * e.desc.max_hrtf_voices);
        e.freeVoices.reset(std::max(e.desc.max_hrtf_voices, 1u));
        for (std::uint32_t i = 0; i < e.desc.max_hrtf_voices; ++i) {
            HrtfVoice& v = e.voices[i];
            v = HrtfVoice{};
            float* cursor = e.voiceMemory.data() + perVoice * i;
            const auto take = [&](std::size_t count) {
                float* p = cursor;
                cursor += count;
                return p;
            };
            v.input = take(2ull * B);
            v.fdlRe = take(P * B);
            v.fdlIm = take(P * B);
            v.filterRe = take(4 * P * B);
            v.filterIm = take(4 * P * B);
            v.delay[0] = take(2ull * kDelayFrames);
            v.delay[1] = take(2ull * kDelayFrames);
            for (float*& band : v.bands) band = take(kReflectionFrames);
            v.splitLow.setCoefficients(designLowpass(400.0f, 0.7071f, e.rate));
            v.splitHigh.setCoefficients(designHighpass(2500.0f, 0.7071f, e.rate));
            e.freeVoices[e.freeVoiceCount++] = e.desc.max_hrtf_voices - 1 - i;
        }
        const std::size_t shChannels = kMaxChannels;
        e.shMix.reset(shChannels * B);
        e.shInput.reset(shChannels * B * 4);
        const std::size_t shPartitions = 2048 / B + 1;
        if (e.C) {
            e.shFdlRe.reset(static_cast<std::size_t>(e.C) * shPartitions * B);
            e.shFdlIm.reset(static_cast<std::size_t>(e.C) * shPartitions * B);
        }
        for (std::uint32_t i = 0; i < std::max(d.max_beds, 1u); ++i) e.beds[i] = RenderBed{};
        e.smoothing = d.smoothing_seconds > 0.0f ? 1.0f - std::exp(-static_cast<float>(B) / (e.rate * d.smoothing_seconds)) : 1.0f;
        fitAirShelves(e);
        ms_head_set_orientation(&e, ms_quat{0, 0, 0, 1});
        *engine = owned.release();
        return MS_OK;
    } catch (const std::bad_alloc&) {
        return MS_ERR_OUT_OF_MEMORY;
    }
}

void ms_engine_destroy(ms_engine* engine) {
    if (!engine) return;
    for (const auto& held : engine->heldSets) ms_hrtf_release(asHandle(held.set));
    delete engine;
}

ms_status ms_engine_commit(ms_engine* engine) {
    if (!engine) return MS_ERR_INVALID_ARGUMENT;
    ms_engine& e = *engine;
    Guard guard(engine);
    collectRetired(e);
    if (!e.currentSet) {
        std::unique_ptr<HrtfSet> fallback = mspatial::hrtf::buildSphereSet(0.0875f, e.rate);
        if (!fallback) return MS_ERR_OUT_OF_MEMORY;
        const ms_status status = selectHrtf(e, fallback.get());
        if (status != MS_OK) return status;
        ms_hrtf_release(asHandle(fallback.release()));
    }
    const std::uint32_t w = e.triple.writeIndex();
    e.snapListener[w] = e.workListener;
    std::copy(e.workSources.begin(), e.workSources.begin() + static_cast<std::ptrdiff_t>(e.highestSource), e.snapSources[w].data());
    std::copy(e.workBeds.begin(), e.workBeds.end(), e.snapBeds[w].data());
    std::copy(e.workBuses.begin(), e.workBuses.end(), e.snapBuses[w].data());
    e.triple.publish();
    // One commit lands in one block: a stereo pair or a layered hit must not start a block apart. Commands that find too
    // little room all wait for the next commit, and only a batch larger than the whole ring is ever split.
    const std::uint32_t capacity = e.commands.capacity();
    while (e.pendingHead < e.pending.size()) {
        const auto count = static_cast<std::uint32_t>(std::min<std::size_t>(e.pending.size() - e.pendingHead, capacity));
        if (!e.commands.pushAll(e.pending.data() + e.pendingHead, count)) break;
        e.pendingHead += count;
    }
    if (e.pendingHead == e.pending.size()) {
        e.pending.clear();
        e.pendingHead = 0;
    }
    return MS_OK;
}

ms_status ms_engine_set_block_callback(ms_engine* engine, ms_block_callback callback, void* user) {
    if (!engine) return MS_ERR_INVALID_ARGUMENT;
    Guard guard(engine);
    queue(*engine, CommandType::BlockCallback, 0, 0, 0, 0, 0.0, reinterpret_cast<void*>(callback), user);
    return MS_OK;
}

uint32_t ms_engine_direct_latency_frames(const ms_engine*) { return 4; }
float ms_engine_sample_rate(const ms_engine* engine) { return engine->rate; }
uint32_t ms_engine_block_frames(const ms_engine* engine) { return engine->B; }

void ms_hrtf_import_options_default(ms_hrtf_import_options* o) {
    *o = ms_hrtf_import_options{};
    o->diffuse_field_eq = 1;
}

ms_status ms_hrtf_load(const char* path, float sampleRate, ms_hrtf** hrtf) {
    if (!path || !hrtf) return MS_ERR_INVALID_ARGUMENT;
    try {
        ms_status status;
        std::unique_ptr<HrtfSet> set = mspatial::hrtf::loadSet(path, status);
        if (!set) return status;
        if (sampleRate > 0.0f && sampleRate != set->sampleRate) set = mspatial::hrtf::resampleSet(*set, sampleRate);
        *hrtf = asHandle(set.release());
        return MS_OK;
    } catch (const std::bad_alloc&) {
        return MS_ERR_OUT_OF_MEMORY;
    }
}

ms_status ms_hrtf_import_sofa(const char* path, const ms_hrtf_import_options* options, ms_hrtf** hrtf) {
    if (!path || !hrtf) return MS_ERR_INVALID_ARGUMENT;
    ms_hrtf_import_options defaults;
    ms_hrtf_import_options_default(&defaults);
    const ms_hrtf_import_options& o = options ? *options : defaults;
    if (o.taps > kMaxFilterTaps) return MS_ERR_INVALID_ARGUMENT;
    try {
        ms_status status;
        const std::unique_ptr<mspatial::hrtf::RawHrirSet> raw = mspatial::hrtf::readSofa(path, status);
        if (!raw) return status;
        mspatial::hrtf::BuildOptions build;
        build.sampleRate = o.sample_rate;
        build.taps = o.taps;
        build.diffuseFieldEq = o.diffuse_field_eq != 0;
        build.headRadius = o.head_radius_m;
        std::unique_ptr<HrtfSet> set = mspatial::hrtf::buildSet(*raw, build);
        if (!set) return MS_ERR_BAD_FILE;
        *hrtf = asHandle(set.release());
        return MS_OK;
    } catch (const std::bad_alloc&) {
        return MS_ERR_OUT_OF_MEMORY;
    }
}

ms_status ms_hrtf_create_sphere(float radius, float sampleRate, ms_hrtf** hrtf) {
    if (!hrtf || !(radius > 0.03f && radius < 0.2f) || !(sampleRate >= 8000.0f && sampleRate <= 192000.0f)) return MS_ERR_INVALID_ARGUMENT;
    try {
        *hrtf = asHandle(mspatial::hrtf::buildSphereSet(radius, sampleRate).release());
        return *hrtf ? MS_OK : MS_ERR_OUT_OF_MEMORY;
    } catch (const std::bad_alloc&) {
        return MS_ERR_OUT_OF_MEMORY;
    }
}

ms_status ms_hrtf_save(const ms_hrtf* hrtf, const char* path) {
    if (!hrtf || !path) return MS_ERR_INVALID_ARGUMENT;
    return mspatial::hrtf::saveSet(*asSet(hrtf), path);
}

void ms_hrtf_get_info(const ms_hrtf* hrtf, ms_hrtf_info* info) {
    const HrtfSet& s = *asSet(hrtf);
    *info = ms_hrtf_info{};
    std::snprintf(info->name, sizeof info->name, "%s", s.name.c_str());
    std::snprintf(info->licence, sizeof info->licence, "%s", s.licence.c_str());
    info->sample_rate = s.sampleRate;
    info->head_radius_m = s.headRadius;
    info->measured_distance_m = s.distance;
    info->directions = s.directions;
    info->taps = s.taps;
}

void ms_hrtf_retain(ms_hrtf* hrtf) { asSet(hrtf)->references.fetch_add(1, std::memory_order_relaxed); }
void ms_hrtf_release(ms_hrtf* hrtf) {
    if (hrtf && asSet(hrtf)->references.fetch_sub(1, std::memory_order_acq_rel) == 1) delete asSet(hrtf);
}

ms_status ms_engine_set_hrtf(ms_engine* engine, ms_hrtf* hrtf) {
    if (!engine || !hrtf) return MS_ERR_INVALID_ARGUMENT;
    Guard guard(engine);
    try {
        return selectHrtf(*engine, asSet(hrtf));
    } catch (const std::bad_alloc&) {
        return MS_ERR_OUT_OF_MEMORY;
    }
}

ms_status ms_listener_set_pose(ms_engine* engine, ms_vec3 position, ms_vec3 velocity, ms_vec3 forward, ms_vec3 up) {
    if (!engine || !finite(position) || !finite(velocity) || !finite(forward) || !finite(up)) return MS_ERR_INVALID_ARGUMENT;
    Guard guard(engine);
    engine->workListener = {position, velocity, forward, up};
    return MS_OK;
}

void ms_head_set_orientation(ms_engine* engine, ms_quat q) {
    if (!engine || !std::isfinite(q.x + q.y + q.z + q.w)) return;
    const Quat n = normalized(Quat(q));
    const auto pack = [](float v) { return static_cast<std::uint64_t>(static_cast<std::uint16_t>(static_cast<std::int16_t>(std::lround(std::clamp(v, -1.0f, 1.0f) * 32767.0f)))); };
    engine->headPacked.store(pack(n.x) | pack(n.y) << 16 | pack(n.z) << 32 | pack(n.w) << 48, std::memory_order_release);
}

ms_status ms_buffer_create(ms_engine* engine, const float* interleaved, uint64_t frames, uint32_t channels, float sampleRate, ms_buffer* buffer) {
    if (!engine || !interleaved || frames == 0 || channels == 0 || channels > kMaxChannels || !buffer || !(sampleRate >= 1000.0f && sampleRate <= 768000.0f))
        return MS_ERR_INVALID_ARGUMENT;
    ms_engine& e = *engine;
    try {
        auto object = std::make_unique<BufferObject>();
        object->channels = channels;
        if (sampleRate == e.rate) {
            object->frames = frames;
            object->data.reset(frames * channels);
            std::memcpy(object->data.data(), interleaved, frames * channels * sizeof(float));
        } else {
            object->frames = resampledLength(frames, sampleRate, e.rate);
            object->data.reset(object->frames * channels);
            std::vector<float> plane(object->frames);
            for (std::uint32_t c = 0; c < channels; ++c) {
                resampleSinc(interleaved + c, frames, channels, sampleRate, e.rate, plane.data(), object->frames);
                for (std::uint64_t i = 0; i < object->frames; ++i) object->data[i * channels + c] = plane[i];
            }
        }
        Guard guard(engine);
        const std::uint32_t id = allocateSlot(e.bufferSlots);
        if (!id) return MS_ERR_LIMIT_REACHED;
        e.buffers[(id & kIndexMask) - 1] = std::move(object);
        *buffer = id;
        return MS_OK;
    } catch (const std::bad_alloc&) {
        return MS_ERR_OUT_OF_MEMORY;
    }
}

ms_status ms_buffer_load_file(ms_engine* engine, const char* path, ms_buffer* buffer) {
    if (!engine || !path || !buffer) return MS_ERR_INVALID_ARGUMENT;
    CFURLRef url = CFURLCreateFromFileSystemRepresentation(nullptr, reinterpret_cast<const UInt8*>(path), static_cast<CFIndex>(std::strlen(path)), false);
    if (!url) return MS_ERR_INVALID_ARGUMENT;
    ExtAudioFileRef file = nullptr;
    const OSStatus opened = ExtAudioFileOpenURL(url, &file);
    CFRelease(url);
    if (opened != noErr || !file) return MS_ERR_NOT_FOUND;
    struct Close {
        ExtAudioFileRef f;
        ~Close() { ExtAudioFileDispose(f); }
    } close{file};
    AudioStreamBasicDescription format{};
    UInt32 size = sizeof format;
    if (ExtAudioFileGetProperty(file, kExtAudioFileProperty_FileDataFormat, &size, &format) != noErr) return MS_ERR_BAD_FILE;
    SInt64 frames = 0;
    size = sizeof frames;
    if (ExtAudioFileGetProperty(file, kExtAudioFileProperty_FileLengthFrames, &size, &frames) != noErr || frames <= 0) return MS_ERR_BAD_FILE;
    const std::uint32_t channels = format.mChannelsPerFrame;
    if (channels == 0 || channels > kMaxChannels) return MS_ERR_UNSUPPORTED;
    AudioStreamBasicDescription client{};
    client.mFormatID = kAudioFormatLinearPCM;
    client.mFormatFlags = kAudioFormatFlagIsFloat | kAudioFormatFlagIsPacked;
    client.mSampleRate = format.mSampleRate;
    client.mChannelsPerFrame = channels;
    client.mBitsPerChannel = 32;
    client.mBytesPerFrame = 4 * channels;
    client.mFramesPerPacket = 1;
    client.mBytesPerPacket = client.mBytesPerFrame;
    if (ExtAudioFileSetProperty(file, kExtAudioFileProperty_ClientDataFormat, sizeof client, &client) != noErr) return MS_ERR_BAD_FILE;
    try {
        std::vector<float> data(static_cast<std::size_t>(frames) * channels);
        std::uint64_t filled = 0;
        while (filled < static_cast<std::uint64_t>(frames)) {
            UInt32 want = static_cast<UInt32>(std::min<std::uint64_t>(static_cast<std::uint64_t>(frames) - filled, 1u << 20));
            AudioBufferList list{};
            list.mNumberBuffers = 1;
            list.mBuffers[0].mNumberChannels = channels;
            list.mBuffers[0].mDataByteSize = want * 4 * channels;
            list.mBuffers[0].mData = data.data() + filled * channels;
            if (ExtAudioFileRead(file, &want, &list) != noErr) return MS_ERR_BAD_FILE;
            if (want == 0) break;
            filled += want;
        }
        if (filled == 0) return MS_ERR_BAD_FILE;
        return ms_buffer_create(engine, data.data(), filled, channels, static_cast<float>(format.mSampleRate), buffer);
    } catch (const std::bad_alloc&) {
        return MS_ERR_OUT_OF_MEMORY;
    }
}

ms_status ms_buffer_destroy(ms_engine* engine, ms_buffer buffer) {
    if (!engine) return MS_ERR_INVALID_ARGUMENT;
    Guard guard(engine);
    std::uint32_t index;
    if (!validSlot(engine->bufferSlots, buffer, index)) return MS_ERR_NOT_FOUND;
    // Kept until the render side says it has let go; the slot itself can be reused at once.
    queue(*engine, CommandType::BufferRetire, 0, 0, 0, 0, 0.0, engine->buffers[index].get());
    engine->retiringBuffers.push_back(std::move(engine->buffers[index]));
    freeSlot(engine->bufferSlots, index);
    return MS_OK;
}

ms_status ms_buffer_get_info(ms_engine* engine, ms_buffer buffer, uint64_t* frames, uint32_t* channels) {
    if (!engine) return MS_ERR_INVALID_ARGUMENT;
    Guard guard(engine);
    std::uint32_t index;
    if (!validSlot(engine->bufferSlots, buffer, index)) return MS_ERR_NOT_FOUND;
    if (frames) *frames = engine->buffers[index]->frames;
    if (channels) *channels = engine->buffers[index]->channels;
    return MS_OK;
}

void ms_source_params_default(ms_source_params* p) {
    *p = ms_source_params{};
    p->forward = {0, 0, -1};
    p->gain = 1.0f;
    p->pitch = 1.0f;
    p->mode = MS_SPATIAL_HRTF;
    p->distance_model = MS_DISTANCE_INVERSE;
    p->ref_distance = 1.0f;
    p->max_distance = 1000.0f;
    p->rolloff = 1.0f;
    p->doppler = 1.0f;
    p->air_absorption = 1.0f;
    p->near_field = 1.0f;
    p->directivity_sharpness = 1.0f;
}

ms_status ms_source_create(ms_engine* engine, ms_source* source) {
    if (!engine || !source) return MS_ERR_INVALID_ARGUMENT;
    ms_engine& e = *engine;
    Guard guard(engine);
    const std::uint32_t id = allocateSlot(e.sourceSlots);
    if (!id) return MS_ERR_LIMIT_REACHED;
    const std::uint32_t index = (id & kIndexMask) - 1;
    ms_source_params_default(&e.workSources[index].params);
    e.workSources[index].reflectionCount = 0;
    e.highestSource = std::max(e.highestSource, index + 1);
    queue(e, CommandType::SourceReset, index);
    *source = id;
    return MS_OK;
}

#define MS_SOURCE_INDEX()                                                     \
    if (!engine) return MS_ERR_INVALID_ARGUMENT;                              \
    ms_engine& e = *engine;                                                   \
    Guard guard(engine);                                                      \
    std::uint32_t index;                                                      \
    if (!validSlot(e.sourceSlots, source, index)) return MS_ERR_NOT_FOUND

ms_status ms_source_destroy(ms_engine* engine, ms_source source) {
    MS_SOURCE_INDEX();
    queue(e, CommandType::SourceReset, index);
    freeSlot(e.sourceSlots, index);
    return MS_OK;
}

ms_status ms_source_set_params(ms_engine* engine, ms_source source, const ms_source_params* params) {
    if (!params || !validParams(*params)) return MS_ERR_INVALID_ARGUMENT;
    MS_SOURCE_INDEX();
    std::uint32_t bus;
    if (params->bus && !validSlot(e.busSlots, params->bus, bus)) return MS_ERR_NOT_FOUND;
    for (ms_bus send : params->send_bus)
        if (send && !validSlot(e.busSlots, send, bus)) return MS_ERR_NOT_FOUND;
    e.workSources[index].params = *params;
    return MS_OK;
}

ms_status ms_source_get_params(ms_engine* engine, ms_source source, ms_source_params* params) {
    if (!params) return MS_ERR_INVALID_ARGUMENT;
    MS_SOURCE_INDEX();
    *params = e.workSources[index].params;
    return MS_OK;
}

ms_status ms_source_set_pose(ms_engine* engine, ms_source source, ms_vec3 position, ms_vec3 velocity) {
    if (!finite(position) || !finite(velocity)) return MS_ERR_INVALID_ARGUMENT;
    MS_SOURCE_INDEX();
    e.workSources[index].params.position = position;
    e.workSources[index].params.velocity = velocity;
    return MS_OK;
}

ms_status ms_source_set_buffer(ms_engine* engine, ms_source source, ms_buffer buffer, uint32_t channel, int loop) {
    MS_SOURCE_INDEX();
    std::uint32_t at;
    if (!validSlot(e.bufferSlots, buffer, at)) return MS_ERR_NOT_FOUND;
    if (channel >= e.buffers[at]->channels) return MS_ERR_INVALID_ARGUMENT;
    queue(e, CommandType::SourceBuffer, index, channel, loop != 0, 0, 0.0, e.buffers[at].get());
    return MS_OK;
}

ms_status ms_source_set_callback(ms_engine* engine, ms_source source, ms_source_callback callback, void* user) {
    if (!callback) return MS_ERR_INVALID_ARGUMENT;
    MS_SOURCE_INDEX();
    queue(e, CommandType::SourceCallback, index, 0, 0, 0, 0.0, reinterpret_cast<void*>(callback), user);
    return MS_OK;
}

ms_status ms_source_play(ms_engine* engine, ms_source source, uint64_t atFrame) {
    MS_SOURCE_INDEX();
    queue(e, CommandType::SourcePlay, index, 0, 0, atFrame);
    return MS_OK;
}

ms_status ms_source_stop(ms_engine* engine, ms_source source) {
    MS_SOURCE_INDEX();
    queue(e, CommandType::SourceStop, index);
    return MS_OK;
}

ms_status ms_source_seek(ms_engine* engine, ms_source source, double frame) {
    if (!std::isfinite(frame) || frame < 0.0) return MS_ERR_INVALID_ARGUMENT;
    MS_SOURCE_INDEX();
    queue(e, CommandType::SourceSeek, index, 0, 0, 0, frame);
    return MS_OK;
}

ms_status ms_source_set_reflections(ms_engine* engine, ms_source source, const ms_reflection_tap* taps, uint32_t count) {
    if (count > MS_MAX_REFLECTION_TAPS || (count && !taps)) return MS_ERR_INVALID_ARGUMENT;
    for (uint32_t i = 0; i < count; ++i) {
        const ms_reflection_tap& t = taps[i];
        if (!std::isfinite(t.delay_seconds) || t.delay_seconds < 0.0f || !finite(t.direction) || !std::isfinite(t.gain[0]) || !std::isfinite(t.gain[1]) ||
            !std::isfinite(t.gain[2]))
            return MS_ERR_INVALID_ARGUMENT;
    }
    MS_SOURCE_INDEX();
    std::copy(taps, taps + count, e.workSources[index].reflections);
    e.workSources[index].reflectionCount = count;
    return MS_OK;
}

ms_status ms_source_get_state(ms_engine* engine, ms_source source, ms_source_state* state) {
    if (!state) return MS_ERR_INVALID_ARGUMENT;
    MS_SOURCE_INDEX();
    const SourcePublished& p = e.published[index];
    state->playing = p.playing.load(std::memory_order_relaxed);
    state->direct = p.direct.load(std::memory_order_relaxed);
    state->position_frames = p.position.load(std::memory_order_relaxed);
    state->level = p.level.load(std::memory_order_relaxed);
    return MS_OK;
}

ms_status ms_bed_create(ms_engine* engine, ms_bed* bed) {
    if (!engine || !bed) return MS_ERR_INVALID_ARGUMENT;
    ms_engine& e = *engine;
    Guard guard(engine);
    if (!e.C) return MS_ERR_UNSUPPORTED;
    const std::uint32_t id = allocateSlot(e.bedSlots);
    if (!id) return MS_ERR_LIMIT_REACHED;
    const std::uint32_t index = (id & kIndexMask) - 1;
    e.workBeds[index] = {Quat{}, 1.0f, 0, 0};
    queue(e, CommandType::BedReset, index);
    *bed = id;
    return MS_OK;
}

#define MS_BED_INDEX()                                                        \
    if (!engine) return MS_ERR_INVALID_ARGUMENT;                              \
    ms_engine& e = *engine;                                                   \
    Guard guard(engine);                                                      \
    std::uint32_t index;                                                      \
    if (!validSlot(e.bedSlots, bed, index)) return MS_ERR_NOT_FOUND

ms_status ms_bed_destroy(ms_engine* engine, ms_bed bed) {
    MS_BED_INDEX();
    queue(e, CommandType::BedReset, index);
    freeSlot(e.bedSlots, index);
    return MS_OK;
}

ms_status ms_bed_set_buffer(ms_engine* engine, ms_bed bed, ms_buffer buffer, int loop) {
    MS_BED_INDEX();
    std::uint32_t at;
    if (!validSlot(e.bufferSlots, buffer, at)) return MS_ERR_NOT_FOUND;
    const std::uint32_t channels = e.buffers[at]->channels;
    const auto order = static_cast<std::uint32_t>(std::lround(std::sqrt(static_cast<double>(channels)))) - 1;
    if ((order + 1) * (order + 1) != channels) return MS_ERR_INVALID_ARGUMENT;
    queue(e, CommandType::BedBuffer, index, order, loop != 0, 0, 0.0, e.buffers[at].get());
    return MS_OK;
}

ms_status ms_bed_set_callback(ms_engine* engine, ms_bed bed, ms_bed_callback callback, void* user, uint32_t order) {
    if (!callback || order > kMaxOrder) return MS_ERR_INVALID_ARGUMENT;
    MS_BED_INDEX();
    queue(e, CommandType::BedCallback, index, order, 0, 0, 0.0, reinterpret_cast<void*>(callback), user);
    return MS_OK;
}

ms_status ms_bed_set_gain(ms_engine* engine, ms_bed bed, float gain) {
    if (!std::isfinite(gain)) return MS_ERR_INVALID_ARGUMENT;
    MS_BED_INDEX();
    e.workBeds[index].gain = gain;
    return MS_OK;
}

ms_status ms_bed_set_bus(ms_engine* engine, ms_bed bed, ms_bus bus) {
    MS_BED_INDEX();
    std::uint32_t busIndex = 0;
    if (bus && !validSlot(e.busSlots, bus, busIndex)) return MS_ERR_NOT_FOUND;
    e.workBeds[index].bus = busIndex;
    return MS_OK;
}

ms_status ms_bed_set_orientation(ms_engine* engine, ms_bed bed, ms_quat orientation, int headLocked) {
    if (!std::isfinite(orientation.x) || !std::isfinite(orientation.y) || !std::isfinite(orientation.z) || !std::isfinite(orientation.w))
        return MS_ERR_INVALID_ARGUMENT;
    MS_BED_INDEX();
    e.workBeds[index].orientation = orientation;
    e.workBeds[index].headLocked = headLocked != 0;
    return MS_OK;
}

ms_status ms_bed_play(ms_engine* engine, ms_bed bed, uint64_t atFrame) {
    MS_BED_INDEX();
    queue(e, CommandType::BedPlay, index, 0, 0, atFrame);
    return MS_OK;
}

ms_status ms_bed_stop(ms_engine* engine, ms_bed bed) {
    MS_BED_INDEX();
    queue(e, CommandType::BedStop, index);
    return MS_OK;
}

ms_status ms_bus_create(ms_engine* engine, ms_bus parent, ms_bus* bus) {
    if (!engine || !bus) return MS_ERR_INVALID_ARGUMENT;
    ms_engine& e = *engine;
    Guard guard(engine);
    std::uint32_t parentIndex;
    if (!validSlot(e.busSlots, parent, parentIndex)) return MS_ERR_NOT_FOUND;
    // Allocated upward only, so parents always have lower indices and the mix order is plain descending. Buses live as
    // long as the engine; a game has a fixed tree.
    for (std::uint32_t i = parentIndex + 1; i < e.busSlots.size(); ++i)
        if (!e.busSlots[i].alive) {
            e.busSlots[i].alive = true;
            // The render side reads parent and used; both are set before the bus can appear in any commit.
            e.buses[i].parent = parentIndex;
            std::atomic_thread_fence(std::memory_order_release);
            std::atomic_ref<bool>(e.buses[i].used).store(true, std::memory_order_release);
            *bus = (static_cast<std::uint32_t>(e.busSlots[i].generation) << 16) | (i + 1);
            return MS_OK;
        }
    return MS_ERR_LIMIT_REACHED;
}

#define MS_BUS_INDEX()                                                        \
    if (!engine) return MS_ERR_INVALID_ARGUMENT;                              \
    ms_engine& e = *engine;                                                   \
    Guard guard(engine);                                                      \
    std::uint32_t index;                                                      \
    if (!validSlot(e.busSlots, bus, index)) return MS_ERR_NOT_FOUND

ms_status ms_bus_set_gain(ms_engine* engine, ms_bus bus, float gain) {
    if (!std::isfinite(gain)) return MS_ERR_INVALID_ARGUMENT;
    MS_BUS_INDEX();
    e.workBuses[index].gain = gain;
    return MS_OK;
}

ms_status ms_bus_set_filter(ms_engine* engine, ms_bus bus, float highpass, float lowpass) {
    if (!std::isfinite(highpass) || !std::isfinite(lowpass) || highpass < 0.0f || lowpass < 0.0f) return MS_ERR_INVALID_ARGUMENT;
    MS_BUS_INDEX();
    e.workBuses[index].highpass = highpass;
    e.workBuses[index].lowpass = lowpass;
    return MS_OK;
}

ms_status ms_bus_set_processor(ms_engine* engine, ms_bus bus, ms_bus_processor processor, void* user) {
    MS_BUS_INDEX();
    queue(e, CommandType::BusProcessor, index, 0, 0, 0, 0.0, reinterpret_cast<void*>(processor), user);
    return MS_OK;
}

void ms_engine_get_stats(const ms_engine* engine, ms_engine_stats* out) {
    const AtomicStats& s = engine->stats;
    const auto relaxed = std::memory_order_relaxed;
    *out = ms_engine_stats{};
    out->frames_rendered = s.frames.load(relaxed);
    out->blocks_rendered = s.blocks.load(relaxed);
    out->sources_playing = s.playing.load(relaxed);
    out->voices_direct = s.direct.load(relaxed);
    out->voices_ambisonic = s.ambisonic.load(relaxed);
    out->commands_pending = engine->commands.pending();
    for (int c = 0; c < 2; ++c) {
        out->peak[c] = s.peak[c].load(relaxed);
        out->rms[c] = s.rms[c].load(relaxed);
    }
    out->overs = s.overs.load(relaxed);
    out->faults = s.faults.load(relaxed);
    out->render_seconds_last = s.last.load(relaxed);
    out->render_seconds_max = s.max.load(relaxed);
}

void ms_engine_get_op_counts(const ms_engine* engine, ms_op_counts* out) {
    const AtomicStats& s = engine->stats;
    const auto relaxed = std::memory_order_relaxed;
    out->blocks = s.opBlocks.load(relaxed);
    out->fft_forward = s.forward.load(relaxed);
    out->fft_inverse = s.inverse.load(relaxed);
    out->spectrum_macs = s.spectrum.load(relaxed);
    out->real_macs = s.real.load(relaxed);
    out->filter_updates = s.updates.load(relaxed);
}

}  // extern "C"
