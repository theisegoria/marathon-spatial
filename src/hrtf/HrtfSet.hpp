// SPDX-License-Identifier: Apache-2.0
// A finished HRTF set: minimum-phase filters plus per-ear onset delays, with an O(1) direction lookup.
#pragma once
#include <atomic>
#include <cmath>
#include <cstdint>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

#include "core/Common.hpp"

namespace mspatial::hrtf {

using core::Aligned;
using core::Vec3;

// Head frame: +x right, +y up, +z behind.
struct RawHrirSet {
    std::vector<Vec3> directions;        // Unit vectors from the head centre.
    std::vector<float> samples;          // [direction][ear (0 left, 1 right)][taps]
    std::vector<float> delays;           // [direction][ear], in samples of sampleRate: a stored onset on top of samples. May be empty.
    std::uint32_t taps{};
    float sampleRate{};
    float distance{};
    std::string name, licence;
};

struct BuildOptions {
    float sampleRate{};  // 0 keeps the input's.
    std::uint32_t taps{};
    bool diffuseFieldEq{true};
    float headRadius{};  // 0 estimates.
};

struct NearFieldEntry {
    float dcDb, hfDb, cornerHz;
};
inline constexpr int kNearRhos = 9, kNearAngles = 19;

struct Neighbours {
    std::uint32_t index[4];
    float weight[4];
};

}  // namespace mspatial::hrtf
namespace mspatial::core { struct AmbisonicDecoder; }
namespace mspatial::hrtf {
using core::AmbisonicDecoder;

class HrtfSet {
public:
    std::string name, licence;
    float sampleRate{}, headRadius{}, distance{};
    std::uint32_t directions{}, taps{};
    Aligned<float> dirs;     // [direction][3]
    Aligned<float> delays;   // [direction][ear], in samples, the smallest in the set being 0
    Aligned<float> filters;  // [direction][ear][taps]
    Aligned<float> quadrature;  // [direction]: solid-angle weights summing to 1
    float maxDelay{};

    // Builds the direction grid and solid-angle weights from dirs.
    void prepareLookup();
    // Finishes a set whose arrays are filled: everything prepareLookup does, and the near-field table.
    void finalize();

    // Modified Shepard interpolation over the four nearest measured directions, with radius from the fifth, so the
    // weights fall continuously to zero as a neighbour leaves the set. Render-side safe: no allocation, O(1).
    void query(Vec3 direction, Neighbours& out) const noexcept;

    // Near-field shelf for an ear, relative to the measured distance. angleToEar in radians, rho = distance / head radius.
    NearFieldEntry nearField(float rho, float angleToEar) const noexcept;

    std::atomic<int> references{1};

    std::shared_ptr<const AmbisonicDecoder> decoder(std::uint32_t order, std::uint32_t blockFrames);

private:
    float fifthDistance(Vec3 direction, Neighbours& out) const noexcept;
    static constexpr std::uint32_t kGrid = 128, kCandidates = 16;
    Aligned<std::uint32_t> grid_;  // [grid][grid][kCandidates]
    NearFieldEntry near_[kNearRhos][kNearAngles]{};
    std::mutex decoderMutex_;
    struct CachedDecoder {
        std::uint32_t order, blockFrames;
        std::shared_ptr<const AmbisonicDecoder> decoder;
    };
    std::vector<CachedDecoder> decoders_;
};

std::unique_ptr<HrtfSet> buildSet(const RawHrirSet& raw, const BuildOptions& options);
std::unique_ptr<HrtfSet> buildSphereSet(float headRadius, float sampleRate);
std::unique_ptr<HrtfSet> resampleSet(const HrtfSet& set, float sampleRate);
ms_status saveSet(const HrtfSet& set, const char* path);
std::unique_ptr<HrtfSet> loadSet(const char* path, ms_status& status);
// Whether the second receiver in a SOFA file is the left ear. The responses decide, since the head shadow cannot lie and some
// published files give ReceiverPosition the opposite sign to their channels. metadataSwapped decides only if they cannot:
// a set with no sources off to the side. positions are Cartesian, [measurement][3]; responses are [measurement][2][taps].
bool receiversSwapped(const float* responses, const float* positions, std::uint32_t measurements, std::uint32_t taps, bool metadataSwapped);
std::unique_ptr<RawHrirSet> readSofa(const char* path, ms_status& status);

// The rigid sphere of Duda and Martens (1998): pressure on the surface relative to free field at the head centre.
// angle is between the source direction and the receiver, rho = source distance / radius, ka = 2 pi f a / c.
void sphereTransfer(double rho, double cosAngle, double ka, double& re, double& im) noexcept;

// Woodworth's high-frequency interaural delay for a far source in the horizontal plane, azimuth 0 to pi / 2.
inline double woodworthItdSeconds(double azimuth, double radius, double speed) {
    return radius / speed * (azimuth + std::sin(azimuth));
}

}  // namespace mspatial::hrtf
