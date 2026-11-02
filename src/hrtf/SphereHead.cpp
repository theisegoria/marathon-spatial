// SPDX-License-Identifier: Apache-2.0
// The analytic rigid sphere: a reference for the tests, a data-free fallback set and the source of the near-field model.
#include <algorithm>
#include <cmath>
#include <complex>
#include <vector>

#include "HrtfSet.hpp"
#include "core/Fft.hpp"

namespace mspatial::hrtf {

using cd = std::complex<double>;

// With an e^{-iwt} convention the exterior field is a sum of outgoing spherical Hankel functions, and the Wronskian of j_n and h_n
// reduces the pressure on a rigid sphere, relative to free field at its centre, to
//   H = -(k r / (k a)^2) Sum (2n + 1) h_n(kr) / h'_n(ka) P_n(x) e^{-ikr}
// Writing h_n(z) = e^{iz} q_n(z) and h'_n(z) = e^{iz} d_n(z), with d_0 = -q_1 and d_n = q_{n-1} - (n + 1) q_n / z:
//   H = -(k r / (k a)^2) e^{-ika} Sum (2n + 1) q_n(kr) / d_n(ka) P_n(x)
// The q_n grow quickly with n at low frequency, and scaling both sequences together at each step keeps their ratio and their
// recurrence intact, so nothing overflows. The result is conjugated into the e^{+iwt} convention signal processing uses, where a
// delay is e^{-iwt}: at angle 0 the wave arrives a / c early, and the result carries the matching phase advance.
void sphereTransfer(double rho, double x, double ka, double& outRe, double& outIm) noexcept {
    const double kr = ka * rho;
    const cd i(0.0, 1.0);
    cd qr0 = 1.0 / (i * kr), qr1 = -(kr + i) / (kr * kr);
    cd qa0 = 1.0 / (i * ka), qa1 = -(ka + i) / (ka * ka);
    double p0 = 1.0, p1 = x;
    cd sum = qr0 / (-qa1);
    for (int n = 1; n < 600; ++n) {
        const cd da = qa0 - (n + 1.0) / ka * qa1;
        const cd term = (2.0 * n + 1.0) * p1 * (qr1 / da);
        sum += term;
        const cd qr2 = (2.0 * n + 1.0) / kr * qr1 - qr0;
        const cd qa2 = (2.0 * n + 1.0) / ka * qa1 - qa0;
        const double p2 = ((2.0 * n + 1.0) * x * p1 - n * p0) / (n + 1.0);
        qr0 = qr1; qr1 = qr2; qa0 = qa1; qa1 = qa2; p0 = p1; p1 = p2;
        const double magnitude = std::abs(qa1);
        if (magnitude > 1e100) {
            const double s = 1.0 / magnitude;
            qr0 *= s; qr1 *= s; qa0 *= s; qa1 *= s;
        }
        if (n > ka + 8 && std::abs(term) < 1e-9 * std::abs(sum)) break;
    }
    const cd H = -(rho / ka) * std::exp(-i * ka) * sum;
    outRe = H.real();
    outIm = -H.imag();
}

std::unique_ptr<HrtfSet> buildSphereSet(float headRadius, float sampleRate) {
    constexpr double speed = 343.0, measureDistance = 1.5;
    constexpr std::uint32_t N = 2048, rawTaps = 1024, preDelay = 128;
    const double rho = measureDistance / headRadius;
    // A table in angle saves evaluating the series at each direction: 1 degree steps, then linear in cos is unnecessary at this
    // resolution for a filter this smooth.
    constexpr int kAngles = 181;
    std::vector<cd> table(static_cast<std::size_t>(kAngles) * (N / 2 + 1));
    for (int a = 0; a < kAngles; ++a) {
        const double x = std::cos(a * M_PI / 180.0);
        for (std::uint32_t k = 1; k <= N / 2; ++k) {
            const double f = k * sampleRate / N;
            const double ka = 2.0 * M_PI * f * headRadius / speed;
            double re, im;
            sphereTransfer(rho, x, ka, re, im);
            table[static_cast<std::size_t>(a) * (N / 2 + 1) + k] = cd(re, im);
        }
        table[static_cast<std::size_t>(a) * (N / 2 + 1)] = cd(1.0, 0.0);
    }
    // A Fibonacci sphere: an even spread with no poles and no symmetry to hide bugs, mirrored so left and right match exactly.
    RawHrirSet raw;
    raw.name = "Rigid sphere";
    raw.licence = "Apache-2.0 (analytic)";
    raw.sampleRate = sampleRate;
    raw.distance = static_cast<float>(measureDistance);
    raw.taps = rawTaps;
    constexpr int halfCount = 512;
    const auto add = [&](Vec3 v) { raw.directions.push_back(core::normalized(v)); };
    const double golden = M_PI * (3.0 - std::sqrt(5.0));
    for (int i = 0; i < halfCount; ++i) {
        const double y = 1.0 - 2.0 * (i + 0.5) / halfCount;
        const double radius = std::sqrt(std::max(0.0, 1.0 - y * y));
        const double phi = golden * i;
        const float px = static_cast<float>(std::abs(std::cos(phi)) * radius);
        const float pz = static_cast<float>(std::sin(phi) * radius);
        add({px, static_cast<float>(y), pz});
        add({-px, static_cast<float>(y), pz});
    }
    for (int a = 0; a < 360; a += 15) {  // Exact horizontal plane directions, which the tests look at.
        const float s = std::sin(static_cast<float>(a) * core::kPi / 180.0f), c = std::cos(static_cast<float>(a) * core::kPi / 180.0f);
        add({s, 0.0f, -c});
    }
    const std::size_t count = raw.directions.size();
    raw.samples.resize(count * 2 * rawTaps);
    core::RealFft fft(N);
    std::vector<float> re(N / 2), im(N / 2), time(N);
    for (std::size_t d = 0; d < count; ++d) {
        for (int ear = 0; ear < 2; ++ear) {
            const double cosine = std::clamp(static_cast<double>(raw.directions[d].x) * (ear == 0 ? -1.0 : 1.0), -1.0, 1.0);
            const double degrees = std::acos(cosine) * 180.0 / M_PI;
            const int a0 = std::min(static_cast<int>(degrees), kAngles - 2);
            const double t = degrees - a0;
            for (std::uint32_t k = 0; k <= N / 2; ++k) {
                const cd value = table[static_cast<std::size_t>(a0) * (N / 2 + 1) + k] * (1.0 - t) +
                                 table[static_cast<std::size_t>(a0 + 1) * (N / 2 + 1) + k] * t;
                // Band limit gently at the top, then delay so causal energy fits. The series' e^{+ika} phase on the near
                // side is an advance, which this delay absorbs.
                const double top = k > N * 0.45 ? 0.5 + 0.5 * std::cos(M_PI * (k - N * 0.45) / (N * 0.05)) : 1.0;
                const cd shifted = value * top * std::polar(1.0, -2.0 * M_PI * k * preDelay / N);
                if (k == 0) re[0] = static_cast<float>(shifted.real());
                else if (k == N / 2) im[0] = static_cast<float>(shifted.real());
                else {
                    re[k] = static_cast<float>(shifted.real());
                    im[k] = static_cast<float>(shifted.imag());
                }
            }
            fft.inverse(re.data(), im.data(), time.data());
            const float scale = 1.0f / N;
            float* out = raw.samples.data() + (d * 2 + static_cast<std::size_t>(ear)) * rawTaps;
            for (std::uint32_t i = 0; i < rawTaps; ++i) out[i] = time[i] * scale;
        }
    }
    BuildOptions options;
    options.diffuseFieldEq = false;
    options.headRadius = headRadius;
    return buildSet(raw, options);
}

}  // namespace mspatial::hrtf
