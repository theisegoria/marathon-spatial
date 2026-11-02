// SPDX-License-Identifier: Apache-2.0
// Higher-order ambisonics, ACN channel order with SN3D normalisation.
#pragma once
#include <cstdint>
#include <memory>

#include "Common.hpp"

namespace mspatial::hrtf { class HrtfSet; }

namespace mspatial::core {

inline constexpr std::uint32_t kMaxOrder = MS_MAX_AMBISONIC_ORDER;
inline constexpr std::uint32_t kMaxChannels = (kMaxOrder + 1) * (kMaxOrder + 1);
inline constexpr std::uint32_t kRotationCells = 680;  // Sum of (2l + 1)^2 up to order 7.

// Ambisonic axes: +x front, +y left, +z up. Head axes: +x right, +y up, +z behind.
inline Vec3 headToAmbisonic(Vec3 v) { return {-v.z, -v.x, v.y}; }

// Writes (order + 1)^2 gains for a unit direction in ambisonic axes. No allocation.
void encodeDirection(Vec3 direction, std::uint32_t order, float* gains) noexcept;

// Spherical-harmonic rotation (Ivanic and Ruedenberg 1996, with their 1998 erratum). The output is one square matrix per
// degree, packed by degree: applying it to the channels of a sound field gives the field rotated by the same rotation.
class ShRotation {
public:
    ShRotation();
    // m is a 3 x 3 rotation in ambisonic axes, row major.
    void compute(const float m[9], std::uint32_t order, float* out) const noexcept;

private:
    struct Coefficients {
        double u, v, w;
    };
    // [l][a + l][b + l]
    Aligned<Coefficients> table_;
    static std::size_t index(int l, int a, int b) noexcept {
        return static_cast<std::size_t>(l) * (2 * kMaxOrder + 1) * (2 * kMaxOrder + 1) +
               static_cast<std::size_t>(a + static_cast<int>(kMaxOrder)) * (2 * kMaxOrder + 1) + static_cast<std::size_t>(b + static_cast<int>(kMaxOrder));
    }
};

// Binaural decode: for each ear and ambisonic channel, a partitioned filter. Built by magnitude least squares (Schoerkhuber,
// Zaunschirm and Hoeldrich 2018): exact least squares where a truncated order can hold the phase, and where it cannot, the
// magnitude is matched and the phase is left free, which keeps timbre and interaural level at low orders.
struct AmbisonicDecoder {
    std::uint32_t order{}, channels{}, blockFrames{}, partitions{};
    std::uint32_t latencyFrames{};
    Aligned<float> re, im;  // [ear][channel][partition][blockFrames]
    const float* spectrumRe(std::uint32_t ear, std::uint32_t channel, std::uint32_t partition) const noexcept {
        return re.data() + (((static_cast<std::size_t>(ear) * channels + channel) * partitions + partition) * blockFrames);
    }
    const float* spectrumIm(std::uint32_t ear, std::uint32_t channel, std::uint32_t partition) const noexcept {
        return im.data() + (((static_cast<std::size_t>(ear) * channels + channel) * partitions + partition) * blockFrames);
    }
};

std::unique_ptr<AmbisonicDecoder> buildDecoder(const hrtf::HrtfSet& set, std::uint32_t order, std::uint32_t blockFrames);

}  // namespace mspatial::core
