// SPDX-License-Identifier: Apache-2.0
// Loudspeaker layouts: panning a direction onto real speakers, and decoding the ambisonic mix to them.
//
// Panning is vector base amplitude panning (Pulkki 1997). A ring of speakers at ear height pans between neighbouring pairs; a layout with
// height pans inside the triangles of the convex hull of the speaker directions. Gaps the speakers leave, usually below the listener and
// often straight above, get imaginary speakers that close the hull; whatever they would have played is dropped and the rest keeps the
// energy. The ambisonic mix is decoded by all-round ambisonic decoding (AllRAD, Zotter and Frank 2012): the field is sampled on a dense,
// even set of virtual speakers with max-rE weights, and each virtual speaker is panned onto the real ones. Nothing here allocates once built.
#pragma once
#include <cstdint>
#include <vector>

#include "Ambisonics.hpp"
#include "Common.hpp"

namespace mspatial::core {

inline constexpr std::uint32_t kMaxSpeakers = MS_MAX_SPEAKERS;

// Head axes, as the rest of the engine: +x right, +y up, +z behind. Azimuth is anticlockwise from the front, so +30 is front left.
Vec3 speakerDirection(float azimuthDegrees, float elevationDegrees);

class SpeakerLayout {
public:
    // False when the layout cannot pan: no speaker that is not LFE, more than kMaxSpeakers, or angles that are not finite.
    bool build(const ms_speaker_layout& layout, std::uint32_t order);
    std::uint32_t count() const { return count_; }
    bool lfe(std::uint32_t speaker) const { return lfe_[speaker]; }

    // Energy-normalised gains on every speaker for a direction in head axes. LFE speakers get none.
    void pan(Vec3 direction, float* gains) const noexcept;
    // Row of the decoder for one speaker: its feed is the sum of these gains times the (order + 1)^2 ambisonic channels.
    const float* decodeRow(std::uint32_t speaker) const { return decode_.data() + static_cast<std::size_t>(speaker) * channels_; }
    std::uint32_t decodeChannels() const { return channels_; }
    // How the stereo master reaches each speaker: its left and right as phantom sources at +30 and -30 degrees.
    const float* stereoIn(std::uint32_t side) const { return stereoIn_[side]; }
    // Each speaker's share of a stereo fold-down: a pan law across +-30 degrees, the rear folding onto its own side.
    const float* stereoOut(std::uint32_t side) const { return stereoOut_[side]; }

private:
    struct Face {
        std::uint32_t a, b, c;
        float inverse[9];  // Rows of the inverse of the matrix whose columns are the three directions.
    };
    struct Pair {
        std::uint32_t a, b;
        float inverse[4];
        float from, arc;  // Azimuth of a and the arc to b, anticlockwise, in radians.
        bool reflex;      // The arc is half a turn or more: no amplitude pair can span it, so it pans by angle.
    };
    enum class Mode { Single, Ring, Sphere };

    void panPoints(Vec3 direction, float* pointGains) const noexcept;

    Mode mode_ = Mode::Single;
    std::uint32_t count_ = 0, channels_ = 0, points_ = 0;
    bool lfe_[kMaxSpeakers]{};
    // Points are the panning speakers followed by any imaginary ones; pointSpeaker maps back, kMaxSpeakers for an imaginary one.
    Vec3 point_[kMaxSpeakers + 2];
    std::uint32_t pointSpeaker_[kMaxSpeakers + 2]{};
    std::vector<Face> faces_;
    std::vector<Pair> pairs_;
    std::vector<float> decode_;
    float stereoIn_[2][kMaxSpeakers]{}, stereoOut_[2][kMaxSpeakers]{};
};

}  // namespace mspatial::core
