#pragma once
// Recorded layers: one-shots, seamless loops and granular clouds over any loaded
// sample, with any parameter (speed, force, a user control) mapped onto gain,
// pitch, filter and grain controls. A LayerInstrument is a BlockSource, so it
// feeds a spatial source as any other mono stream, and it can drive a modal voice
// through ContactKind::Source.
//
// Construction allocates. render, setParameter, trigger and release do not, and
// the output is deterministic for any render block size. One thread renders, and
// setParameter, trigger and release may come from another; they apply at the next
// control tick, every 32 frames.
#include "marathon_spatial/physics/AudioAsset.hpp"
#include "marathon_spatial/physics/Source.hpp"
#include <cstdint>
#include <memory>
#include <vector>

namespace mspatial::physics {
using ParameterId=std::uint32_t;
inline constexpr ParameterId MaximumParameters=64;
inline constexpr ParameterId ParameterSpeed=0,ParameterForce=1,ParameterFirstUser=16;

// inputMinimum..inputMaximum to outputMinimum..outputMaximum. The normalized input
// is clamped, raised to exponent, then mixed linearly, or geometrically when
// logarithmic is set (for hertz, or for gain in constant dB per unit input).
struct ParameterCurve {
    double inputMinimum{},inputMaximum{1},outputMinimum{},outputMaximum{1};
    double exponent{1};
    bool logarithmic{};
};
[[nodiscard]] double evaluateCurve(const ParameterCurve& curve,double input) noexcept;

// Gain and Pitch multiply the layer's base value and each other. The rest replace it.
enum class LayerTarget : std::uint8_t {Gain,Pitch,LowPassHz,GrainDensity,GrainSeconds,GrainPosition,GrainPositionJitter,GrainPitchJitter};
struct ParameterMapping {
    ParameterId parameter{};
    LayerTarget target{LayerTarget::Gain};
    ParameterCurve curve{};
    double smoothingSeconds{.03};
};
enum class LayerKind : std::uint8_t {OneShot,Loop,Granular};
struct SampleLayer {
    LayerKind kind{LayerKind::OneShot};
    std::shared_ptr<const MonoSample> sample;
    double gain{1},pitch{1};
    double lowPassHz{};            // 0: no filter
    std::vector<ParameterMapping> mappings;
    // Loop: the sample's loop region, or the whole sample. The tail is
    // equal-power crossfaded into the head.
    double crossfadeSeconds{.05};
    // Granular
    double grainDensityPerSecond{40},grainSeconds{.05};
    double grainPosition{.5},grainPositionJitter{.2}; // fractions of the sample
    double grainPitchJitter{.02};  // fraction of the pitch, uniform +-
};
struct LayerInstrumentSettings {
    double sampleRateHz{48000};
    std::uint32_t seed{1};
    std::uint32_t maximumGrains{64}; // per granular layer
};
struct LayerStats {std::uint64_t grainsStarted{},grainsRefused{};};
class LayerInstrument final : public BlockSource {
public:
    LayerInstrument(std::vector<SampleLayer> layers,const LayerInstrumentSettings& settings);
    ~LayerInstrument() override;
    LayerInstrument(LayerInstrument&&) noexcept;
    LayerInstrument& operator=(LayerInstrument&&) noexcept;
    void setParameter(ParameterId parameter,double value) noexcept;
    void trigger() noexcept;                    // starts every layer from its beginning
    void release(double seconds) noexcept;      // fades out, then finishes
    void render(std::span<float> mono) noexcept override;
    [[nodiscard]] bool finished() const noexcept override;
    [[nodiscard]] const LayerStats& stats() const noexcept;
private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};
}
