// SPDX-License-Identifier: Apache-2.0
// Between the engine's speaker channels and a Core Audio device: where the device's channels stand, and writing channels into its buffers.
// Pure functions of their arguments, with no device access, so the core's tests reach them and the output only calls them.
#pragma once
#include <CoreAudioTypes/CoreAudioTypes.h>

#include <cstdint>

#include "marathon_spatial/marathon_spatial.h"

namespace mspatial::core {

// Speaker positions from a device's channel layout, in the device's channel order, for deviceChannels channels. Channels the layout names by
// label take that label's usual position (ITU-R BS.2051 angles); channels given by coordinates take those. A channel with no position (an
// unknown or discrete label, a second channel at a position already taken, or one past the layout) is marked lfe: it pans nothing and stays
// silent. False when the layout cannot be read or places no speaker.
bool speakerLayoutFromChannelLayout(const AudioChannelLayout& layout, std::uint32_t deviceChannels, ms_speaker_layout& out);

// Writes channels[c] (frames samples each) into device channel c of output, whether the device's buffers are interleaved or one per channel.
// Device channels past count are silenced. Real-time safe.
void writeDeviceChannels(const float* const* channels, std::uint32_t count, std::uint32_t frames, AudioBufferList& output) noexcept;

}  // namespace mspatial::core
