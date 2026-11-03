// An impulse response from an energy field: first-order ambisonics (ACN W, Y, Z, X; SN3D; world frame mapped as the engine maps AmbiX:
// front -Z, left -X, up +Y).
//
// Each of the field's twelve arrival directions gets its own noise, shaped band by band to the energy the field holds for that direction
// and moment, then encoded into first order. Independent noise per direction is what keeps a diffuse tail diffuse: W, X, Y and Z come out
// mutually uncorrelated in a diffuse field, while a strong direction still pulls the response toward it. The noise is fixed once, so a
// response rebuilt from a slightly different field differs only in its envelope and swaps in without a change of texture.
//
// The twelve noises are one band-limited noise per band read at twelve circular offsets far longer than any band's correlation time, so
// they are uncorrelated in practice and the bank stays at one noise per band.
#pragma once

#include "acoustics/EnergyField.hpp"

#include <vector>

namespace mspatial::acoustics {

inline constexpr std::uint32_t kAmbisonicChannels = 4;

// Where the response starts and stops, in seconds of field time.
struct ReverbIrShape {
    float fadeInStart = 0.f, fadeInEnd = 0.f;    // Rises from zero: the early part belongs to the exact image-source taps.
    float fadeOutStart = 0.f, fadeOutEnd = 0.f;  // Falls to zero: the feedback network carries the tail from here.
    std::uint32_t advance = 0;                   // Samples dropped from the front, to make up for the convolver's latency.
};

class ReverbIrSynth {
public:
    // Allocates the noise bank. maxSamples bounds every response.
    void prepare(float sampleRate, std::uint32_t maxSamples);
    std::uint32_t maxSamples() const { return length_; }
    float sampleRate() const { return sampleRate_; }

    // Writes kAmbisonicChannels channels of maxSamples() each and returns the samples that matter (the rest are zero). The field is
    // smoothed first: more in time as it gets later, since a decaying field changes slowly and its late bins hold the fewest rays.
    std::uint32_t build(const EnergyField& field, const ReverbIrShape& shape, float* const* out);

    // The smoothed field build() used last, per direction, band and bin. For tests and level matching.
    const EnergyField& smoothed() const { return smoothed_; }

private:
    void smooth(const EnergyField& field);

    float sampleRate_ = 48000.f;
    std::uint32_t length_ = 0;
    std::vector<float> noise_;      // [band][length_], unit mean square.
    std::vector<float> envelope_;   // One channel of envelope.
    std::vector<float> window_, position_, amplitude_;
    std::vector<float> signal_;     // One direction's signal.
    EnergyField smoothed_;
    std::vector<float> totals_;
};

}  // namespace mspatial::acoustics
