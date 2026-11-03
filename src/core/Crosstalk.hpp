// SPDX-License-Identifier: Apache-2.0
// Crosstalk cancellation: binaural sound on two loudspeakers.
//
// A binaural signal is meant for the ears, and headphones deliver it there. Loudspeakers do not: each ear also hears the far speaker, a
// little later and duller round the head, and that crosstalk collapses the image to the line between the speakers. A canceller feeds each
// speaker filtered copies of both channels so that, at the ears, the crosstalk cancels (Atal and Schroeder, 1960s; Cooper and Bauck's
// transaural stereo).
//
// The plant is a spherical head between two speakers: Brown and Duda's head-shadow filter for each speaker and ear, Woodworth's delay round
// the head, and the distance from each speaker. The head may sit off centre and turn, as a head tracker reports; the plant is then no
// longer symmetric, so the canceller is a full two-by-two inverse. It is normalised by the near ear's response at the centre, so it removes
// crosstalk and leaves the speakers' own sound alone.
//
// The design follows published least-squares methods:
//  - Tikhonov regularisation (Kirkeby and Nelson's fast deconvolution, 1998) caps the boost where the inverse would blow up. Its weight
//    depends on frequency, as in their frequency-dependent regularisation: looser in the middle of the band, tighter towards its edges,
//    where cancellation holds least against a moving head and a boost colours the most.
//  - The inverse is fitted over several head positions around the given one at once, a least-squares sum over their plants. That gives
//    up a few decibels at the exact point for a much wider sweet spot.
//  - Only crosstalk inside a band, by default 200 Hz to 7 kHz, is cancelled; outside it the canceller is plain stereo by construction.
//  - A gentle correction, at most 6 dB and smoothed over a third of an octave, flattens what the canceller leaves on the near ear.
// The filters add 5.3 ms.
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
    std::uint32_t latencyFrames = 0;   // The modelling delay that keeps the inverse causal.
    std::vector<float> impulse[2][2];  // [speaker][channel]: what each input channel feeds into each speaker.
    Aligned<float> re, im;             // [speaker][channel][partition][blockFrames], packed as the engine's block FFT gives them.
};

// The plant at a frequency for the head where the parameters put it: [ear][speaker], relative to the near ear at the centre.
void crosstalkPlant(const ms_crosstalk_params& params, double hz, double speedOfSound, std::complex<double> out[2][2]);

// Designs the canceller for a rate and the engine's block size. nullptr for parameters out of range.
std::unique_ptr<CrosstalkFilter> designCrosstalk(const ms_crosstalk_params& params, float rate, std::uint32_t blockFrames, double speedOfSound);

}  // namespace mspatial::core
