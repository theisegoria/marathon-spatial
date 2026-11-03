// SPDX-License-Identifier: Apache-2.0
// Crosstalk cancellation: binaural sound on two loudspeakers.
//
// A binaural signal is meant for the ears, and headphones deliver it there. Loudspeakers do not: each ear also hears the far speaker, a
// little later and duller round the head, and that crosstalk collapses the image to the line between the speakers. A canceller feeds each
// speaker a filtered, inverted copy of the other channel that meets the crosstalk at the far ear and removes it (Atal and Schroeder, 1960s;
// Cooper and Bauck's transaural stereo).
//
// The plant is a spherical head between two speakers at +-span: Brown and Duda's head-shadow filter for each ear, with Woodworth's
// interaural delay. It is normalised to the near ear, so only crosstalk is removed and the speakers' own sound is left alone. A symmetric
// plant splits into the sum and difference of the two channels (the shuffler), each inverted on its own with Tikhonov regularisation
// (Kirkeby and Nelson's fast deconvolution), which caps the boost where the inverse would blow up: in the bass, where the two paths differ
// least, and wherever the delayed crosstalk turns in phase. Only crosstalk inside a band, by default 200 Hz to 7 kHz, is cancelled; outside
// it the canceller is plain stereo, because there the cancellation cannot hold against a head that moves. The filters add 5.3 ms.
#pragma once
#include <complex>
#include <cstdint>
#include <memory>
#include <vector>

#include "Common.hpp"

namespace mspatial::core {

inline constexpr std::uint32_t kMaxCrosstalkTaps = 4096;

struct CrosstalkFilter {
    std::uint32_t taps = 0, partitions = 0, blockFrames = 0;
    std::uint32_t latencyFrames = 0;  // The modelling delay that keeps the inverse causal.
    std::vector<float> sum, difference;  // Impulse responses for the sum and difference channels.
    Aligned<float> re, im;               // [channel][partition][blockFrames], packed as the engine's block FFT gives them.
};

// The plant: what reaches the far ear from a speaker, relative to the near ear, at a frequency. Includes the interaural delay.
std::complex<double> crosstalkPlant(const ms_crosstalk_params& params, double hz, double speedOfSound);

// Designs the canceller for a rate and the engine's block size. nullptr for parameters out of range.
std::unique_ptr<CrosstalkFilter> designCrosstalk(const ms_crosstalk_params& params, float rate, std::uint32_t blockFrames, double speedOfSound);

}  // namespace mspatial::core
