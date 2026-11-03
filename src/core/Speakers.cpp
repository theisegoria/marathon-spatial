// SPDX-License-Identifier: Apache-2.0
#include "Speakers.hpp"

#include <algorithm>
#include <array>
#include <cmath>

namespace mspatial::core {
namespace {

constexpr float kRingElevation = 10.0f;      // Degrees. A layout whose speakers all sit within this of ear height pans as a ring.
constexpr float kHullEpsilon = 1e-4f;
constexpr std::uint32_t kDecoderPoints = 2030;  // Virtual speakers of the AllRAD decoder: far more than order 7 needs.

float legendre(std::uint32_t n, float x) {
    if (n == 0) return 1.0f;
    float previous = 1.0f, current = x;
    for (std::uint32_t k = 2; k <= n; ++k) {
        const float next = ((2.0f * k - 1.0f) * x * current - (k - 1.0f) * previous) / static_cast<float>(k);
        previous = current;
        current = next;
    }
    return current;
}

Vec3 fibonacci(std::uint32_t index, std::uint32_t count) {
    const float z = 1.0f - (2.0f * static_cast<float>(index) + 1.0f) / static_cast<float>(count);
    const float r = std::sqrt(std::max(0.0f, 1.0f - z * z));
    const float phi = 2.0f * kPi * std::fmod(static_cast<float>(index) * 0.61803398875f, 1.0f);
    return {r * std::cos(phi), z, r * std::sin(phi)};
}

bool invert3(Vec3 a, Vec3 b, Vec3 c, float out[9]) {
    // Columns a, b, c. The inverse's rows are the cross products over the determinant.
    const Vec3 r0 = cross(b, c), r1 = cross(c, a), r2 = cross(a, b);
    const float det = dot(a, r0);
    if (std::fabs(det) < 1e-6f) return false;
    const float k = 1.0f / det;
    const Vec3 rows[3] = {r0, r1, r2};
    for (int i = 0; i < 3; ++i) {
        out[3 * i] = rows[i].x * k;
        out[3 * i + 1] = rows[i].y * k;
        out[3 * i + 2] = rows[i].z * k;
    }
    return true;
}

}  // namespace

Vec3 speakerDirection(float azimuthDegrees, float elevationDegrees) {
    const float az = azimuthDegrees * kPi / 180.0f, el = elevationDegrees * kPi / 180.0f;
    return {-std::sin(az) * std::cos(el), std::sin(el), -std::cos(az) * std::cos(el)};
}

bool SpeakerLayout::build(const ms_speaker_layout& layout, std::uint32_t order) {
    if (layout.count == 0 || layout.count > kMaxSpeakers || order > kMaxOrder) return false;
    count_ = layout.count;
    points_ = 0;
    faces_.clear();
    pairs_.clear();
    bool ring = true;
    for (std::uint32_t s = 0; s < count_; ++s) {
        const ms_speaker& speaker = layout.speakers[s];
        if (!std::isfinite(speaker.azimuth_degrees) || !std::isfinite(speaker.elevation_degrees)) return false;
        lfe_[s] = speaker.lfe != 0;
        if (lfe_[s]) continue;
        const float elevation = std::clamp(speaker.elevation_degrees, -90.0f, 90.0f);
        const Vec3 d = speakerDirection(speaker.azimuth_degrees, elevation);
        for (std::uint32_t p = 0; p < points_; ++p)
            if (dot(d, point_[p]) > 0.99985f) return false;  // Two speakers within a degree of each other: nothing can tell them apart.
        point_[points_] = d;
        pointSpeaker_[points_++] = s;
        ring = ring && std::fabs(elevation) <= kRingElevation;
    }
    if (points_ == 0) return false;
    const std::uint32_t real = points_;
    mode_ = real == 1 ? Mode::Single : ring ? Mode::Ring : Mode::Sphere;

    if (mode_ == Mode::Ring) {
        // Neighbours by azimuth, all the way round.
        std::array<std::uint32_t, kMaxSpeakers> sorted;
        std::array<float, kMaxSpeakers> azimuth;
        for (std::uint32_t p = 0; p < real; ++p) {
            sorted[p] = p;
            azimuth[p] = std::atan2(-point_[p].x, -point_[p].z);
        }
        std::sort(sorted.begin(), sorted.begin() + real, [&](std::uint32_t a, std::uint32_t b) { return azimuth[a] < azimuth[b]; });
        for (std::uint32_t i = 0; i < real; ++i) {
            const std::uint32_t a = sorted[i], b = sorted[(i + 1) % real];
            float arc = azimuth[b] - azimuth[a];
            if (arc <= 0.0f) arc += 2.0f * kPi;
            const float ax = point_[a].x, az = point_[a].z, bx = point_[b].x, bz = point_[b].z;
            const float det = ax * bz - bx * az;
            Pair pair{a, b, {}, azimuth[a], arc, arc >= kPi - 1e-3f || std::fabs(det) < 1e-6f};
            if (!pair.reflex) {
                pair.inverse[0] = bz / det;
                pair.inverse[1] = -bx / det;
                pair.inverse[2] = -az / det;
                pair.inverse[3] = ax / det;
            }
            pairs_.push_back(pair);
        }
    } else if (mode_ == Mode::Sphere) {
        // Imaginary speakers close the hull where the layout leaves a pole open.
        bool below = false, above = false;
        for (std::uint32_t p = 0; p < real; ++p) {
            below = below || point_[p].y < -0.5f;
            above = above || point_[p].y > 0.87f;
        }
        if (!below) {
            point_[points_] = {0, -1, 0};
            pointSpeaker_[points_++] = kMaxSpeakers;
        }
        if (!above) {
            point_[points_] = {0, 1, 0};
            pointSpeaker_[points_++] = kMaxSpeakers;
        }
        // Convex hull by brute force, which a few dozen points allow. A face may be a polygon of coplanar speakers (symmetric layouts
        // make trapezoids): it is found once and fanned into triangles.
        std::vector<std::vector<std::uint32_t>> seen;
        for (std::uint32_t i = 0; i < points_; ++i)
            for (std::uint32_t j = i + 1; j < points_; ++j)
                for (std::uint32_t k = j + 1; k < points_; ++k) {
                    Vec3 n = cross(point_[j] - point_[i], point_[k] - point_[i]);
                    if (length(n) < 1e-5f) continue;
                    n = normalized(n);
                    const float offset = dot(n, point_[i]);
                    bool positive = false, negative = false;
                    std::vector<std::uint32_t> plane;
                    for (std::uint32_t m = 0; m < points_; ++m) {
                        const float side = dot(n, point_[m]) - offset;
                        if (side > kHullEpsilon) positive = true;
                        else if (side < -kHullEpsilon) negative = true;
                        else plane.push_back(m);
                    }
                    if (positive && negative) continue;
                    if (std::find(seen.begin(), seen.end(), plane) != seen.end()) continue;
                    seen.push_back(plane);
                    Vec3 centre{};
                    for (std::uint32_t m : plane) centre = centre + point_[m];
                    centre = centre * (1.0f / static_cast<float>(plane.size()));
                    const Vec3 u = normalized(point_[plane[0]] - centre), v = cross(n, u);
                    std::sort(plane.begin(), plane.end(), [&](std::uint32_t a, std::uint32_t b) {
                        const Vec3 da = point_[a] - centre, db = point_[b] - centre;
                        return std::atan2(dot(da, v), dot(da, u)) < std::atan2(dot(db, v), dot(db, u));
                    });
                    for (std::size_t t = 1; t + 1 < plane.size(); ++t) {
                        Face face{plane[0], plane[t], plane[t + 1], {}};
                        if (invert3(point_[face.a], point_[face.b], point_[face.c], face.inverse)) faces_.push_back(face);
                    }
                }
        if (faces_.empty()) return false;
    }

    // The stereo master as phantom sources at +-30 degrees; one speaker alone takes both at half power.
    float gains[kMaxSpeakers];
    for (std::uint32_t side = 0; side < 2; ++side) {
        pan(speakerDirection(side == 0 ? 30.0f : -30.0f, 0.0f), gains);
        for (std::uint32_t s = 0; s < count_; ++s) stereoIn_[side][s] = mode_ == Mode::Single ? gains[s] * std::sqrt(0.5f) : gains[s];
    }
    for (std::uint32_t s = 0; s < count_; ++s) {
        stereoOut_[0][s] = stereoOut_[1][s] = 0.0f;
        if (lfe_[s]) continue;
        const Vec3 d = speakerDirection(layout.speakers[s].azimuth_degrees, 0.0f);
        const float azimuth = std::atan2(-d.x, -d.z) * 180.0f / kPi;
        const float position = std::clamp(azimuth / 30.0f, -1.0f, 1.0f);
        const float angle = (1.0f - position) * 0.25f * kPi;
        stereoOut_[0][s] = std::cos(angle);
        stereoOut_[1][s] = std::sin(angle);
    }

    // AllRAD: sample the field with max-rE weights on even virtual speakers, pan each onto the real ones, then scale so that a plane wave
    // carries, on average over directions, the energy of a direct source.
    channels_ = (order + 1) * (order + 1);
    decode_.assign(static_cast<std::size_t>(count_) * channels_, 0.0f);
    const float spread = std::cos(137.9f * kPi / 180.0f / (static_cast<float>(order) + 1.51f));
    std::vector<float> weight(channels_);
    for (std::uint32_t n = 0; n <= order; ++n)
        for (std::uint32_t c = n * n; c < (n + 1) * (n + 1); ++c) weight[c] = (2.0f * n + 1.0f) * legendre(n, spread) / kDecoderPoints;
    float encode[kMaxChannels];
    for (std::uint32_t k = 0; k < kDecoderPoints; ++k) {
        const Vec3 d = fibonacci(k, kDecoderPoints);
        pan(d, gains);
        encodeDirection(headToAmbisonic(d), order, encode);
        for (std::uint32_t s = 0; s < count_; ++s) {
            if (gains[s] == 0.0f) continue;
            float* row = decode_.data() + static_cast<std::size_t>(s) * channels_;
            for (std::uint32_t c = 0; c < channels_; ++c) row[c] += gains[s] * weight[c] * encode[c];
        }
    }
    double energy = 0.0;
    for (std::uint32_t k = 0; k < kDecoderPoints; ++k) {
        encodeDirection(headToAmbisonic(fibonacci(k, kDecoderPoints)), order, encode);
        for (std::uint32_t s = 0; s < count_; ++s) {
            const float* row = decodeRow(s);
            float feed = 0.0f;
            for (std::uint32_t c = 0; c < channels_; ++c) feed += row[c] * encode[c];
            energy += static_cast<double>(feed) * feed;
        }
    }
    const float scale = energy > 0.0 ? static_cast<float>(1.0 / std::sqrt(energy / kDecoderPoints)) : 0.0f;
    for (float& value : decode_) value *= scale;
    return true;
}

void SpeakerLayout::panPoints(Vec3 direction, float* pointGains) const noexcept {
    for (std::uint32_t p = 0; p < points_; ++p) pointGains[p] = 0.0f;
    if (mode_ == Mode::Single) {
        pointGains[0] = 1.0f;
        return;
    }
    if (mode_ == Mode::Ring) {
        const float x = direction.x, z = direction.z;
        if (x * x + z * z < 1e-8f) {  // Straight up or down: every speaker of the ring alike.
            for (std::uint32_t p = 0; p < points_; ++p) pointGains[p] = 1.0f;
            return;
        }
        // The pair whose arc holds the direction. An arc of half a turn or more (behind a front-only layout) pans by angle instead.
        const float angle = std::atan2(-x, -z);
        for (const Pair& pair : pairs_) {
            float into = angle - pair.from;
            while (into < 0.0f) into += 2.0f * kPi;
            while (into >= 2.0f * kPi) into -= 2.0f * kPi;
            if (into > pair.arc) continue;
            if (pair.reflex) {
                const float t = into / pair.arc;
                pointGains[pair.a] = std::cos(0.5f * kPi * t);
                pointGains[pair.b] = std::sin(0.5f * kPi * t);
            } else {
                pointGains[pair.a] = std::max(pair.inverse[0] * x + pair.inverse[1] * z, 0.0f);
                pointGains[pair.b] = std::max(pair.inverse[2] * x + pair.inverse[3] * z, 0.0f);
            }
            return;
        }
        for (std::uint32_t p = 0; p < points_; ++p) pointGains[p] = 1.0f;
        return;
    }
    float best = -INFINITY, g[3]{};
    const Face* chosen = nullptr;
    for (const Face& face : faces_) {
        const float a = face.inverse[0] * direction.x + face.inverse[1] * direction.y + face.inverse[2] * direction.z;
        const float b = face.inverse[3] * direction.x + face.inverse[4] * direction.y + face.inverse[5] * direction.z;
        const float c = face.inverse[6] * direction.x + face.inverse[7] * direction.y + face.inverse[8] * direction.z;
        const float worst = std::min(a, std::min(b, c));
        if (worst > best) {
            best = worst;
            chosen = &face;
            g[0] = a;
            g[1] = b;
            g[2] = c;
        }
    }
    if (!chosen) return;
    pointGains[chosen->a] = std::max(g[0], 0.0f);
    pointGains[chosen->b] = std::max(g[1], 0.0f);
    pointGains[chosen->c] = std::max(g[2], 0.0f);
}

void SpeakerLayout::pan(Vec3 direction, float* gains) const noexcept {
    for (std::uint32_t s = 0; s < count_; ++s) gains[s] = 0.0f;
    if (points_ == 0) return;
    float pointGains[kMaxSpeakers + 2];
    panPoints(normalized(direction), pointGains);
    float energy = 0.0f;
    for (std::uint32_t p = 0; p < points_; ++p) {
        if (pointSpeaker_[p] == kMaxSpeakers) continue;  // Imaginary: what it would play is dropped.
        gains[pointSpeaker_[p]] = pointGains[p];
        energy += pointGains[p] * pointGains[p];
    }
    if (energy < 1e-12f) {
        // Everything fell on an imaginary speaker, straight at an open pole: spread it over every real one.
        for (std::uint32_t p = 0; p < points_; ++p)
            if (pointSpeaker_[p] != kMaxSpeakers) gains[pointSpeaker_[p]] = 1.0f, energy += 1.0f;
    }
    const float scale = 1.0f / std::sqrt(energy);
    for (std::uint32_t s = 0; s < count_; ++s) gains[s] *= scale;
}

}  // namespace mspatial::core
