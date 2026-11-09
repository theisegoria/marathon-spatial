// SPDX-License-Identifier: Apache-2.0
#include "DeviceChannels.hpp"

#include "Speakers.hpp"

#include <AudioToolbox/AudioFormat.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <vector>

namespace mspatial::core {
namespace {

struct Position {
    float azimuth, elevation;  // Degrees, azimuth anticlockwise from the front.
    bool lfe;
};

bool positionOf(const AudioChannelDescription& d, Position& out) {
    const auto at = [&](float azimuth, float elevation) {
        out = {azimuth, elevation, false};
        return true;
    };
    switch (d.mChannelLabel) {
    case kAudioChannelLabel_Left:
    case kAudioChannelLabel_LeftTotal:
    case kAudioChannelLabel_HeadphonesLeft:
    case kAudioChannelLabel_BinauralLeft: return at(30, 0);
    case kAudioChannelLabel_Right:
    case kAudioChannelLabel_RightTotal:
    case kAudioChannelLabel_HeadphonesRight:
    case kAudioChannelLabel_BinauralRight: return at(-30, 0);
    case kAudioChannelLabel_Center:
    case kAudioChannelLabel_Mono: return at(0, 0);
    case kAudioChannelLabel_LFEScreen:
    case kAudioChannelLabel_LFE2: out = {0, 0, true}; return true;
    case kAudioChannelLabel_LeftSurround: return at(110, 0);
    case kAudioChannelLabel_RightSurround: return at(-110, 0);
    case kAudioChannelLabel_LeftCenter: return at(15, 0);
    case kAudioChannelLabel_RightCenter: return at(-15, 0);
    case kAudioChannelLabel_CenterSurround:
    case kAudioChannelLabel_CenterSurroundDirect: return at(180, 0);
    case kAudioChannelLabel_LeftSurroundDirect: return at(90, 0);
    case kAudioChannelLabel_RightSurroundDirect: return at(-90, 0);
    case kAudioChannelLabel_RearSurroundLeft: return at(150, 0);
    case kAudioChannelLabel_RearSurroundRight: return at(-150, 0);
    case kAudioChannelLabel_LeftWide: return at(60, 0);
    case kAudioChannelLabel_RightWide: return at(-60, 0);
    case kAudioChannelLabel_TopCenterSurround: return at(0, 90);
    case kAudioChannelLabel_VerticalHeightLeft: return at(30, 30);
    case kAudioChannelLabel_VerticalHeightCenter: return at(0, 30);
    case kAudioChannelLabel_VerticalHeightRight: return at(-30, 30);
    case kAudioChannelLabel_TopBackLeft: return at(135, 30);
    case kAudioChannelLabel_TopBackCenter: return at(180, 30);
    case kAudioChannelLabel_TopBackRight: return at(-135, 30);
    case kAudioChannelLabel_LeftTopMiddle: return at(90, 45);
    case kAudioChannelLabel_RightTopMiddle: return at(-90, 45);
    case kAudioChannelLabel_LeftTopRear: return at(135, 45);
    case kAudioChannelLabel_CenterTopRear: return at(180, 45);
    case kAudioChannelLabel_RightTopRear: return at(-135, 45);
    case kAudioChannelLabel_UseCoordinates:
        if (d.mChannelFlags & kAudioChannelFlags_SphericalCoordinates)  // Core Audio's azimuth is positive to the right.
            return at(-d.mCoordinates[kAudioChannelCoordinates_Azimuth], d.mCoordinates[kAudioChannelCoordinates_Elevation]);
        if (d.mChannelFlags & kAudioChannelFlags_RectangularCoordinates) {
            const float x = d.mCoordinates[kAudioChannelCoordinates_LeftRight], y = d.mCoordinates[kAudioChannelCoordinates_BackFront];
            const float z = d.mCoordinates[kAudioChannelCoordinates_DownUp];
            if (x == 0.0f && y == 0.0f && z == 0.0f) return false;
            constexpr float kDegrees = 57.29577951f;
            return at(std::atan2(-x, y) * kDegrees, std::atan2(z, std::hypot(x, y)) * kDegrees);
        }
        return false;
    default: return false;
    }
}

// A layout named by a tag or a bitmap, expanded into descriptions; the layout itself when it has them. Null when it cannot be expanded.
const AudioChannelLayout* describe(const AudioChannelLayout& layout, std::vector<unsigned char>& expanded) {
    if (layout.mChannelLayoutTag == kAudioChannelLayoutTag_UseChannelDescriptions) return &layout;
    const bool bitmap = layout.mChannelLayoutTag == kAudioChannelLayoutTag_UseChannelBitmap;
    const AudioFormatPropertyID property = bitmap ? kAudioFormatProperty_ChannelLayoutForBitmap : kAudioFormatProperty_ChannelLayoutForTag;
    const void* specifier = bitmap ? static_cast<const void*>(&layout.mChannelBitmap) : static_cast<const void*>(&layout.mChannelLayoutTag);
    const UInt32 specifierSize = bitmap ? sizeof layout.mChannelBitmap : sizeof layout.mChannelLayoutTag;
    UInt32 size = 0;
    if (AudioFormatGetPropertyInfo(property, specifierSize, specifier, &size) != noErr || size < offsetof(AudioChannelLayout, mChannelDescriptions)) return nullptr;
    expanded.resize(size);
    if (AudioFormatGetProperty(property, specifierSize, specifier, &size, expanded.data()) != noErr) return nullptr;
    return reinterpret_cast<const AudioChannelLayout*>(expanded.data());
}

}  // namespace

bool stereoSubFromChannelLayout(const AudioChannelLayout* layout, std::uint32_t deviceChannels, std::uint32_t& left, std::uint32_t& right,
                                std::uint32_t& sub) {
    constexpr std::uint32_t kNone = ~0u;
    left = right = sub = kNone;
    std::vector<unsigned char> expanded;
    const AudioChannelLayout* described = layout ? describe(*layout, expanded) : nullptr;
    if (described) {
        const std::uint32_t n = std::min<std::uint32_t>(deviceChannels, described->mNumberChannelDescriptions);
        for (std::uint32_t c = 0; c < n; ++c) {
            Position p{};
            if (!positionOf(described->mChannelDescriptions[c], p)) continue;
            if (p.lfe) {
                if (sub == kNone) sub = c;
            } else if (p.elevation == 0.0f && std::fabs(p.azimuth - 30.0f) < 1.0f) {
                if (left == kNone) left = c;
            } else if (p.elevation == 0.0f && std::fabs(p.azimuth + 30.0f) < 1.0f) {
                if (right == kNone) right = c;
            }
        }
        if (left != kNone && right != kNone && sub != kNone) return true;
        if (left != kNone || right != kNone || sub != kNone) return false;  // Labelled, but not as a 2.1 set.
    }
    if (deviceChannels == 3) {  // Nothing labelled: three channels are taken as left, right, LFE.
        left = 0;
        right = 1;
        sub = 2;
        return true;
    }
    left = right = sub = kNone;
    return false;
}

bool speakerLayoutFromChannelLayout(const AudioChannelLayout& layout, std::uint32_t deviceChannels, ms_speaker_layout& out) {
    out = ms_speaker_layout{};
    std::vector<unsigned char> expanded;
    const AudioChannelLayout* described = describe(layout, expanded);
    if (!described) return false;
    out.count = std::min<std::uint32_t>(deviceChannels, MS_MAX_SPEAKERS);
    bool any = false;
    for (std::uint32_t c = 0; c < out.count; ++c) {
        ms_speaker& speaker = out.speakers[c];
        speaker = {0, 0, 1};  // No position until one is found.
        Position p{};
        if (c >= described->mNumberChannelDescriptions || !positionOf(described->mChannelDescriptions[c], p) || p.lfe) continue;
        if (!std::isfinite(p.azimuth) || !std::isfinite(p.elevation)) continue;
        // Two channels within a degree of each other cannot be told apart; the first keeps the place, as the layout builder requires.
        const float elevation = std::clamp(p.elevation, -90.0f, 90.0f);
        const Vec3 here = speakerDirection(p.azimuth, elevation);
        bool taken = false;
        for (std::uint32_t k = 0; k < c && !taken; ++k)
            taken = !out.speakers[k].lfe &&
                    dot(here, speakerDirection(out.speakers[k].azimuth_degrees, out.speakers[k].elevation_degrees)) > 0.99985f;
        if (taken) continue;
        speaker = {p.azimuth, elevation, 0};
        any = true;
    }
    return any;
}

void writeDeviceChannels(const float* const* channels, std::uint32_t count, std::uint32_t frames, AudioBufferList& output) noexcept {
    std::uint32_t first = 0;  // Device channel index of this buffer's first channel.
    for (UInt32 b = 0; b < output.mNumberBuffers; ++b) {
        AudioBuffer& buffer = output.mBuffers[b];
        const UInt32 stride = std::max<UInt32>(buffer.mNumberChannels, 1);
        auto* data = static_cast<float*>(buffer.mData);
        const UInt32 length = std::min<UInt32>(frames, buffer.mDataByteSize / (sizeof(float) * stride));
        if (!data) {
            first += stride;
            continue;
        }
        if (stride == 1) {
            if (first < count) std::memcpy(data, channels[first], length * sizeof(float));
            else std::memset(data, 0, buffer.mDataByteSize);
        } else {
            for (UInt32 k = 0; k < stride; ++k) {
                const std::uint32_t channel = first + k;
                if (channel < count)
                    for (UInt32 t = 0; t < length; ++t) data[t * stride + k] = channels[channel][t];
                else
                    for (UInt32 t = 0; t < length; ++t) data[t * stride + k] = 0.0f;
            }
        }
        first += stride;
    }
}

}  // namespace mspatial::core
