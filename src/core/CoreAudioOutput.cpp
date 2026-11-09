// SPDX-License-Identifier: Apache-2.0
// Live output on a Core Audio device. The HAL's I/O thread calls render directly: time-constraint priority, one block per
// cycle, no intermediate thread or queue. The thread joins the device's audio workgroup, so the scheduler sees its deadline
// and keeps it on a performance core. Nothing here runs until ms_output_start.
#include <CoreAudio/CoreAudio.h>
#include <os/workgroup.h>
#include <unistd.h>

#include <algorithm>
#include <atomic>
#include <cerrno>
#include <memory>
#include <vector>

#include "DeviceChannels.hpp"
#include "WorkgroupJoin.hpp"
#include "marathon_spatial/marathon_spatial_output.h"

struct ms_output {
    ms_engine* engine{};
    AudioObjectID device{kAudioObjectUnknown};
    AudioDeviceIOProcID proc{};
    os_workgroup_t workgroup{};
    os_workgroup_join_token_s token{};
    mspatial::core::WorkgroupJoin join;  // Belongs to the device thread; reset only while that thread is stopped.
    std::vector<float> left, right;
    std::vector<float> channelMemory;  // One maxFrames run per engine channel, for a speaker engine on a device with room for it.
    std::vector<float*> channels;
    std::uint32_t deviceChannels{}, engineChannels{2};
    // 2.1: a stereo engine's left, right and subwoofer, mapped onto the device's channels (the rest point at silence).
    std::int32_t subChannel{-1};
    std::vector<float> sub, zeros;
    std::vector<const float*> deviceMap;
    bool running{};
    bool listening{};
    std::uint32_t maxFrames{};
    std::atomic<bool> joined{false}, wantLeave{false}, alive{true}, rateChanged{false};
    std::atomic<std::uint64_t> callbacks{0}, overruns{0};
    std::atomic<std::uint32_t> bufferFrames{0}, joinedOnce{0}, latency{0};
};

namespace {

AudioObjectPropertyAddress address(AudioObjectPropertySelector selector, AudioObjectPropertyScope scope = kAudioObjectPropertyScopeGlobal) {
    return {selector, scope, kAudioObjectPropertyElementMain};
}

template <class T>
OSStatus getProperty(AudioObjectID object, AudioObjectPropertySelector selector, AudioObjectPropertyScope scope, T& value) {
    const AudioObjectPropertyAddress a = address(selector, scope);
    UInt32 size = sizeof value;
    return AudioObjectGetPropertyData(object, &a, 0, nullptr, &size, &value);
}

OSStatus deviceListener(AudioObjectID, UInt32 count, const AudioObjectPropertyAddress* addresses, void* user) {
    auto* o = static_cast<ms_output*>(user);
    for (UInt32 i = 0; i < count; ++i) {
        if (addresses[i].mSelector == kAudioDeviceProcessorOverload) o->overruns.fetch_add(1, std::memory_order_relaxed);
        if (addresses[i].mSelector == kAudioDevicePropertyDeviceIsAlive) o->alive.store(false, std::memory_order_relaxed);
        if (addresses[i].mSelector == kAudioDevicePropertyNominalSampleRate) o->rateChanged.store(true, std::memory_order_relaxed);
    }
    return noErr;
}

constexpr AudioObjectPropertySelector kWatched[] = {kAudioDeviceProcessorOverload, kAudioDevicePropertyDeviceIsAlive, kAudioDevicePropertyNominalSampleRate};

OSStatus ioProc(AudioObjectID, const AudioTimeStamp*, const AudioBufferList*, const AudioTimeStamp*, AudioBufferList* output, const AudioTimeStamp*, void* user) {
    auto* o = static_cast<ms_output*>(user);
    if (!o->joined.load(std::memory_order_relaxed) && !o->wantLeave.load(std::memory_order_relaxed) && o->workgroup) {
        const mspatial::core::WorkgroupJoin::Result result = o->join.step([o] { return os_workgroup_join(o->workgroup, &o->token); });
        if (result == mspatial::core::WorkgroupJoin::Result::Joined) {
            o->joined.store(true, std::memory_order_release);
            o->joinedOnce.store(1, std::memory_order_relaxed);
        } else if (result == mspatial::core::WorkgroupJoin::Result::AlreadyMember) {  // The HAL has already put this thread in the group.
            o->joinedOnce.store(1, std::memory_order_relaxed);
        }
    }
    if (o->wantLeave.load(std::memory_order_acquire) && o->joined.load(std::memory_order_relaxed)) {
        os_workgroup_leave(o->workgroup, &o->token);  // A thread must leave from itself.
        o->joined.store(false, std::memory_order_release);
    }
    o->callbacks.fetch_add(1, std::memory_order_relaxed);
    if (!output || output->mNumberBuffers == 0) return noErr;
    const AudioBuffer& first = output->mBuffers[0];
    const UInt32 stride = std::max<UInt32>(first.mNumberChannels, 1);
    const UInt32 frames = first.mDataByteSize / (sizeof(float) * stride);
    o->bufferFrames.store(frames, std::memory_order_relaxed);
    if (frames > o->maxFrames || o->rateChanged.load(std::memory_order_relaxed)) {  // Silence is the only honest answer.
        for (UInt32 b = 0; b < output->mNumberBuffers; ++b) std::memset(output->mBuffers[b].mData, 0, output->mBuffers[b].mDataByteSize);
        return noErr;
    }
    if (o->subChannel < 0 && o->engineChannels > 2) {  // A channel per speaker, each to its own device channel.
        ms_engine_render_speakers(o->engine, o->channels.data(), frames);
        mspatial::core::writeDeviceChannels(o->channels.data(), o->engineChannels, frames, *output);
        return noErr;
    }
    float* left = o->left.data();
    float* right = o->right.data();
    if (o->subChannel >= 0) {
        ms_engine_render_2_1(o->engine, left, right, o->sub.data(), frames);
        mspatial::core::writeDeviceChannels(o->deviceMap.data(), static_cast<std::uint32_t>(o->deviceMap.size()), frames, *output);
        return noErr;
    }
    ms_engine_render(o->engine, left, right, frames);
    if (stride >= 2) {
        auto* out = static_cast<float*>(first.mData);
        std::memset(out, 0, first.mDataByteSize);
        for (UInt32 i = 0; i < frames; ++i) {
            out[i * stride] = left[i];
            out[i * stride + 1] = right[i];
        }
        for (UInt32 b = 1; b < output->mNumberBuffers; ++b) std::memset(output->mBuffers[b].mData, 0, output->mBuffers[b].mDataByteSize);
    } else {
        for (UInt32 b = 0; b < output->mNumberBuffers; ++b) {
            if (b < 2) std::memcpy(output->mBuffers[b].mData, b == 0 ? left : right, frames * sizeof(float));
            else std::memset(output->mBuffers[b].mData, 0, output->mBuffers[b].mDataByteSize);
        }
        if (output->mNumberBuffers == 1) {  // A mono device gets the mid signal.
            auto* out = static_cast<float*>(first.mData);
            for (UInt32 i = 0; i < frames; ++i) out[i] = 0.5f * (left[i] + right[i]);
        }
    }
    return noErr;
}

std::uint32_t outputChannels(AudioObjectID device) {
    const AudioObjectPropertyAddress a = address(kAudioDevicePropertyStreamConfiguration, kAudioObjectPropertyScopeOutput);
    UInt32 size = 0;
    if (AudioObjectGetPropertyDataSize(device, &a, 0, nullptr, &size) != noErr || size < sizeof(AudioBufferList)) return 0;
    std::vector<unsigned char> bytes(size);
    if (AudioObjectGetPropertyData(device, &a, 0, nullptr, &size, bytes.data()) != noErr) return 0;
    const auto* list = reinterpret_cast<const AudioBufferList*>(bytes.data());
    std::uint32_t total = 0;
    for (UInt32 b = 0; b < list->mNumberBuffers; ++b) total += list->mBuffers[b].mNumberChannels;
    return total;
}

ms_status resolveDevice(const char* uid, AudioObjectID& device) {
    if (!uid) {
        const AudioObjectPropertyAddress a = address(kAudioHardwarePropertyDefaultOutputDevice);
        UInt32 size = sizeof device;
        if (AudioObjectGetPropertyData(kAudioObjectSystemObject, &a, 0, nullptr, &size, &device) != noErr || device == kAudioObjectUnknown) return MS_ERR_DEVICE;
        return MS_OK;
    }
    CFStringRef string = CFStringCreateWithCString(nullptr, uid, kCFStringEncodingUTF8);
    if (!string) return MS_ERR_INVALID_ARGUMENT;
    const AudioObjectPropertyAddress a = address(kAudioHardwarePropertyTranslateUIDToDevice);
    UInt32 size = sizeof device;
    const OSStatus status = AudioObjectGetPropertyData(kAudioObjectSystemObject, &a, sizeof string, &string, &size, &device);
    CFRelease(string);
    return status == noErr && device != kAudioObjectUnknown ? MS_OK : MS_ERR_NOT_FOUND;
}

}  // namespace

extern "C" {

ms_status ms_output_default_sample_rate(float* sampleRate) {
    if (!sampleRate) return MS_ERR_INVALID_ARGUMENT;
    AudioObjectID device = kAudioObjectUnknown;
    if (const ms_status status = resolveDevice(nullptr, device); status != MS_OK) return status;
    Float64 rate = 0.0;
    if (getProperty(device, kAudioDevicePropertyNominalSampleRate, kAudioObjectPropertyScopeGlobal, rate) != noErr) return MS_ERR_DEVICE;
    *sampleRate = static_cast<float>(rate);
    return MS_OK;
}

namespace {
// The device's preferred output channel layout, or empty.
std::vector<unsigned char> preferredLayout(AudioObjectID device) {
    const AudioObjectPropertyAddress a = address(kAudioDevicePropertyPreferredChannelLayout, kAudioObjectPropertyScopeOutput);
    UInt32 size = 0;
    if (!AudioObjectHasProperty(device, &a) || AudioObjectGetPropertyDataSize(device, &a, 0, nullptr, &size) != noErr ||
        size < offsetof(AudioChannelLayout, mChannelDescriptions))
        return {};
    std::vector<unsigned char> bytes(size);
    if (AudioObjectGetPropertyData(device, &a, 0, nullptr, &size, bytes.data()) != noErr) return {};
    return bytes;
}
// left, right, sub as device channels, or false.
bool stereoSub(AudioObjectID device, std::uint32_t channels, std::uint32_t& left, std::uint32_t& right, std::uint32_t& sub) {
    if (channels < 3 || channels > 64) return false;
    const auto bytes = preferredLayout(device);
    return mspatial::core::stereoSubFromChannelLayout(bytes.empty() ? nullptr : reinterpret_cast<const AudioChannelLayout*>(bytes.data()), channels,
                                                      left, right, sub);
}
}  // namespace

ms_status ms_output_device_sub_channel(const char* uid, int32_t* channel) {
    if (!channel) return MS_ERR_INVALID_ARGUMENT;
    *channel = -1;
    AudioObjectID device = kAudioObjectUnknown;
    if (const ms_status status = resolveDevice(uid, device); status != MS_OK) return status;
    std::uint32_t left, right, sub;
    if (stereoSub(device, outputChannels(device), left, right, sub)) *channel = static_cast<int32_t>(sub);
    return MS_OK;
}

ms_status ms_output_device_layout(const char* uid, ms_speaker_layout* layout) {
    if (!layout) return MS_ERR_INVALID_ARGUMENT;
    AudioObjectID device = kAudioObjectUnknown;
    if (const ms_status status = resolveDevice(uid, device); status != MS_OK) return status;
    const std::uint32_t channels = outputChannels(device);
    if (channels == 0) return MS_ERR_DEVICE;
    const AudioObjectPropertyAddress a = address(kAudioDevicePropertyPreferredChannelLayout, kAudioObjectPropertyScopeOutput);
    UInt32 size = 0;
    if (AudioObjectHasProperty(device, &a) && AudioObjectGetPropertyDataSize(device, &a, 0, nullptr, &size) == noErr &&
        size >= offsetof(AudioChannelLayout, mChannelDescriptions)) {
        std::vector<unsigned char> bytes(size);
        if (AudioObjectGetPropertyData(device, &a, 0, nullptr, &size, bytes.data()) == noErr &&
            mspatial::core::speakerLayoutFromChannelLayout(*reinterpret_cast<const AudioChannelLayout*>(bytes.data()), channels, *layout))
            return MS_OK;
    }
    if (channels == 2) return ms_speaker_layout_preset(MS_SPEAKERS_STEREO, layout);
    return MS_ERR_NOT_FOUND;
}

ms_status ms_output_create(ms_engine* engine, const char* uid, ms_output** output) {
    if (!engine || !output) return MS_ERR_INVALID_ARGUMENT;
    auto o = std::make_unique<ms_output>();
    o->engine = engine;
    if (const ms_status status = resolveDevice(uid, o->device); status != MS_OK) return status;

    Float64 rate = 0.0;
    if (getProperty(o->device, kAudioDevicePropertyNominalSampleRate, kAudioObjectPropertyScopeGlobal, rate) != noErr) return MS_ERR_DEVICE;
    if (static_cast<float>(rate) != ms_engine_sample_rate(engine)) return MS_ERR_RATE_MISMATCH;

    // A hardware cycle of one block if the device allows it, so each callback renders exactly one block.
    const UInt32 block = ms_engine_block_frames(engine);
    AudioValueRange range{};
    if (getProperty(o->device, kAudioDevicePropertyBufferFrameSizeRange, kAudioObjectPropertyScopeGlobal, range) != noErr) return MS_ERR_DEVICE;
    UInt32 frames = static_cast<UInt32>(std::clamp<Float64>(block, range.mMinimum, range.mMaximum));
    const AudioObjectPropertyAddress bufferSize = address(kAudioDevicePropertyBufferFrameSize);
    AudioObjectSetPropertyData(o->device, &bufferSize, 0, nullptr, sizeof frames, &frames);  // A refusal only costs a little FIFO.

    o->maxFrames = static_cast<std::uint32_t>(std::min<Float64>(range.mMaximum, 16384.0));
    o->left.resize(o->maxFrames);
    o->right.resize(o->maxFrames);

    UInt32 deviceLatency = 0, safety = 0;
    (void)getProperty(o->device, kAudioDevicePropertyLatency, kAudioObjectPropertyScopeOutput, deviceLatency);
    (void)getProperty(o->device, kAudioDevicePropertySafetyOffset, kAudioObjectPropertyScopeOutput, safety);
    UInt32 streamsBytes = 0;
    const AudioObjectPropertyAddress streams = address(kAudioDevicePropertyStreams, kAudioObjectPropertyScopeOutput);
    UInt32 streamLatency = 0;
    if (AudioObjectGetPropertyDataSize(o->device, &streams, 0, nullptr, &streamsBytes) == noErr && streamsBytes >= sizeof(AudioStreamID)) {
        std::vector<AudioStreamID> ids(streamsBytes / sizeof(AudioStreamID));
        if (AudioObjectGetPropertyData(o->device, &streams, 0, nullptr, &streamsBytes, ids.data()) == noErr) {
            AudioStreamBasicDescription format{};
            if (getProperty(ids[0], kAudioStreamPropertyVirtualFormat, kAudioObjectPropertyScopeGlobal, format) != noErr ||
                format.mBitsPerChannel != 32 || !(format.mFormatFlags & kAudioFormatFlagIsFloat))
                return MS_ERR_UNSUPPORTED;
            (void)getProperty(ids[0], kAudioStreamPropertyLatency, kAudioObjectPropertyScopeGlobal, streamLatency);
        }
    } else {
        return MS_ERR_DEVICE;  // No output streams.
    }
    o->latency.store(deviceLatency + safety + streamLatency, std::memory_order_relaxed);
    o->deviceChannels = outputChannels(o->device);
    // A speaker engine gets a channel per speaker when the device has them all; anything else plays stereo, folded down if need be.
    const std::uint32_t speakers = ms_engine_output_channels(engine);
    if (speakers > 2 && o->deviceChannels >= speakers) {
        o->engineChannels = speakers;
        o->channelMemory.resize(static_cast<std::size_t>(speakers) * o->maxFrames);
        for (std::uint32_t c = 0; c < speakers; ++c) o->channels.push_back(o->channelMemory.data() + static_cast<std::size_t>(c) * o->maxFrames);
    }
    std::uint32_t left, right, sub;
    if (speakers <= 2 && stereoSub(o->device, o->deviceChannels, left, right, sub)) {
        o->subChannel = static_cast<std::int32_t>(sub);
        o->engineChannels = 3;
        o->sub.resize(o->maxFrames);
        o->zeros.assign(o->maxFrames, 0.0f);
        o->deviceMap.assign(o->deviceChannels, o->zeros.data());
        o->deviceMap[left] = o->left.data();
        o->deviceMap[right] = o->right.data();
        o->deviceMap[sub] = o->sub.data();
    }

    os_workgroup_t workgroup = nullptr;
    if (getProperty(o->device, kAudioDevicePropertyIOThreadOSWorkgroup, kAudioObjectPropertyScopeGlobal, workgroup) == noErr) o->workgroup = workgroup;
    if (AudioDeviceCreateIOProcID(o->device, ioProc, o.get(), &o->proc) != noErr) {
        if (o->workgroup) os_release(o->workgroup);
        return MS_ERR_DEVICE;
    }
    *output = o.release();
    return MS_OK;
}

ms_status ms_output_start(ms_output* o) {
    if (!o) return MS_ERR_INVALID_ARGUMENT;
    if (o->running) return MS_OK;
    o->wantLeave.store(false, std::memory_order_release);
    o->join.reset();
    for (AudioObjectPropertySelector selector : kWatched) {
        const AudioObjectPropertyAddress a = address(selector);
        AudioObjectAddPropertyListener(o->device, &a, deviceListener, o);
    }
    o->listening = true;
    if (AudioDeviceStart(o->device, o->proc) != noErr) return MS_ERR_DEVICE;
    o->running = true;
    return MS_OK;
}

ms_status ms_output_stop(ms_output* o) {
    if (!o) return MS_ERR_INVALID_ARGUMENT;
    if (o->running) {
        o->wantLeave.store(true, std::memory_order_release);
        for (int i = 0; i < 40 && o->joined.load(std::memory_order_acquire); ++i) usleep(5000);
        AudioDeviceStop(o->device, o->proc);
        o->running = false;
    }
    if (o->listening) {
        for (AudioObjectPropertySelector selector : kWatched) {
            const AudioObjectPropertyAddress a = address(selector);
            AudioObjectRemovePropertyListener(o->device, &a, deviceListener, o);
        }
        o->listening = false;
    }
    return MS_OK;
}

void ms_output_destroy(ms_output* o) {
    if (!o) return;
    ms_output_stop(o);
    AudioDeviceDestroyIOProcID(o->device, o->proc);
    if (o->workgroup) os_release(o->workgroup);
    delete o;
}

void ms_output_get_stats(const ms_output* o, ms_output_stats* s) {
    *s = ms_output_stats{};
    s->callbacks = o->callbacks.load(std::memory_order_relaxed);
    s->overruns = o->overruns.load(std::memory_order_relaxed);
    s->buffer_frames = o->bufferFrames.load(std::memory_order_relaxed);
    s->workgroup_joined = o->joinedOnce.load(std::memory_order_relaxed);
    s->device_alive = o->alive.load(std::memory_order_relaxed);
    s->latency_frames = o->latency.load(std::memory_order_relaxed);
    s->device_channels = o->deviceChannels;
    s->engine_channels = o->engineChannels;
    s->sub_channel = o->subChannel;
}

}  // extern "C"
