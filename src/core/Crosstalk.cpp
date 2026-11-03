// SPDX-License-Identifier: Apache-2.0
#include "Crosstalk.hpp"

#include <algorithm>
#include <array>
#include <cmath>

#include "Fft.hpp"

namespace mspatial::core {
namespace {

constexpr double kPiD = 3.14159265358979323846;
using Complex = std::complex<double>;
using Matrix = std::array<std::array<Complex, 2>, 2>;

// Brown and Duda's head shadow at the ear for a source theta radians off that ear's axis.
Complex shadow(double omega, double theta, double radius, double c) {
    constexpr double alphaMin = 0.1, thetaMin = 150.0 * kPiD / 180.0;
    const double alpha = (1.0 + alphaMin / 2.0) + (1.0 - alphaMin / 2.0) * std::cos(theta / thetaMin * kPiD);
    const double w0 = c / radius;
    return Complex(1.0, alpha * omega / (2.0 * w0)) / Complex(1.0, omega / (2.0 * w0));
}

// Woodworth's delay at the ear relative to the head's centre: shorter for an ear facing the source, round the sphere for one that is not.
double earDelay(double theta, double radius, double c) {
    return theta < kPiD / 2.0 ? -radius / c * std::cos(theta) : radius / c * (theta - kPiD / 2.0);
}

// 1 inside [low, high], falling to 0 over an octave either side along a raised cosine.
double bandWeight(double hz, double low, double high) {
    const auto edge = [](double octaves) { return octaves >= 1.0 ? 0.0 : 0.5 * (1.0 + std::cos(kPiD * octaves)); };
    if (hz <= 0.0) return 0.0;
    if (hz < low) return edge(std::log2(low / hz));
    if (hz > high) return edge(std::log2(hz / high));
    return 1.0;
}

// 0 in the middle of the band, an octave in from either edge, rising to 1 at the edges and beyond: where the regulariser tightens.
double edgeShape(double hz, double low, double high) {
    if (hz <= low || hz >= high) return 1.0;
    const double octaves = std::min(std::log2(hz / low), std::log2(high / hz));
    return octaves >= 1.0 ? 0.0 : 0.5 * (1.0 + std::cos(kPiD * octaves));
}

struct Head {
    double x, z, yaw;  // Metres right and back of the centre, and radians turned to the left.
};

// The plant before normalisation: [ear][speaker], with the speakers at distance and span and the head where it is.
Matrix rawPlant(const ms_crosstalk_params& p, const Head& head, double hz, double c) {
    const double span = p.speaker_span_degrees * kPiD / 180.0, a = p.head_radius_m, D = p.listener_distance_m, omega = 2.0 * kPiD * hz;
    const double speakerX[2] = {-D * std::sin(span), D * std::sin(span)}, speakerZ = -D * std::cos(span);
    // Facing -Z turned by yaw towards -X; the right ear's axis is the head's right.
    const double rightX = std::cos(head.yaw), rightZ = -std::sin(head.yaw);
    Matrix h{};
    for (int speaker = 0; speaker < 2; ++speaker) {
        const double ux = speakerX[speaker] - head.x, uz = speakerZ - head.z, r = std::hypot(ux, uz);
        for (int ear = 0; ear < 2; ++ear) {
            const double side = ear == 0 ? -1.0 : 1.0;
            const double theta = std::acos(std::clamp(side * (ux * rightX + uz * rightZ) / r, -1.0, 1.0));
            h[ear][speaker] = D / r * shadow(omega, theta, a, c) * std::polar(1.0, -omega * (r / c + earDelay(theta, a, c)));
        }
    }
    return h;
}

Matrix multiply(const Matrix& a, const Matrix& b) {
    Matrix out{};
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j) out[i][j] = a[i][0] * b[0][j] + a[i][1] * b[1][j];
    return out;
}

Matrix adjoint(const Matrix& a) { return {{{std::conj(a[0][0]), std::conj(a[1][0])}, {std::conj(a[0][1]), std::conj(a[1][1])}}}; }

Matrix inverse(const Matrix& a) {
    const Complex det = a[0][0] * a[1][1] - a[0][1] * a[1][0];
    return {{{a[1][1] / det, -a[0][1] / det}, {-a[1][0] / det, a[0][0] / det}}};
}

bool valid(const ms_crosstalk_params& p) {
    if (!(p.speaker_span_degrees >= 2.0f && p.speaker_span_degrees <= 60.0f) || !(p.head_radius_m >= 0.05f && p.head_radius_m <= 0.15f)) return false;
    if (!(p.max_boost_db >= 3.0f && p.max_boost_db <= 30.0f) || !(p.edge_boost_db >= 3.0f && p.edge_boost_db <= p.max_boost_db)) return false;
    if (!(p.low_hz >= 20.0f) || !(p.high_hz > p.low_hz * 2.0f) || !(p.listener_distance_m >= 0.2f && p.listener_distance_m <= 5.0f)) return false;
    if (!(std::fabs(p.head_x_m) <= 0.5f) || !(std::fabs(p.head_z_m) <= 0.5f) || !(std::fabs(p.head_yaw_degrees) <= 60.0f)) return false;
    if (!(p.robust_radius_m >= 0.0f && p.robust_radius_m <= 0.1f)) return false;
    // The head, wherever the design puts it, stays clear of both speakers.
    const double span = p.speaker_span_degrees * kPiD / 180.0, D = p.listener_distance_m;
    for (double sx : {-D * std::sin(span), D * std::sin(span)})
        if (std::hypot(sx - p.head_x_m, -D * std::cos(span) - p.head_z_m) < p.head_radius_m + p.robust_radius_m + 0.05) return false;
    return true;
}

}  // namespace

void crosstalkPlant(const ms_crosstalk_params& p, double hz, double c, Complex out[2][2]) {
    const Matrix h = rawPlant(p, {p.head_x_m, p.head_z_m, p.head_yaw_degrees * kPiD / 180.0}, hz, c);
    const Complex reference = rawPlant(p, {0.0, 0.0, 0.0}, hz, c)[0][0];
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j) out[i][j] = h[i][j] / reference;
}

std::unique_ptr<CrosstalkFilter> designCrosstalk(const ms_crosstalk_params& p, float rate, std::uint32_t B, double c) {
    if (!valid(p) || !(rate > 0.0f) || B == 0) return nullptr;
    auto f = std::make_unique<CrosstalkFilter>();
    f->taps = rate > 64000.0f ? 4096 : 2048;  // About 43 ms: fine enough in the bass that the band edge stays put.
    f->blockFrames = B;
    f->partitions = (f->taps + B - 1) / B;
    // The modelling delay. The regularised inverse rings on both sides of its main tap, longest in the bass, where the two ears' paths differ
    // least; an eighth of the filter, 5.3 ms at 48 kHz, holds that ringing, where a sixteenth gives up about 7 dB of cancellation.
    f->latencyFrames = f->taps / 8;

    // The head positions the inverse is fitted over: the one given, and a cross around it in the horizontal plane.
    const double yaw = p.head_yaw_degrees * kPiD / 180.0, rho = p.robust_radius_m;
    std::vector<Head> heads{{p.head_x_m, p.head_z_m, yaw}};
    if (rho > 0.0)
        for (const auto& [dx, dz] : {std::pair{rho, 0.0}, std::pair{-rho, 0.0}, std::pair{0.0, rho}, std::pair{0.0, -rho}})
            heads.push_back({p.head_x_m + dx, p.head_z_m + dz, yaw});
    const double K = static_cast<double>(heads.size());

    // Design on a grid four times the filter. The band is applied to the plant: inside it the canceller sees the plant, outside it plain
    // stereo, so the canceller is plain stereo there by construction and nothing has to be blended (blending a finished inverse with
    // plain stereo crosses two responses of very different phase at the band edge, which rings far longer than the filter).
    // Each bin's inverse is the regularised least-squares fit over the heads, C = (sum H* H + K beta I)^-1 sum H*. A single plant's
    // singular values s come out as s / (s^2 + beta), at most 1 / (2 sqrt(beta)): beta is set by the boost cap, looser in the middle of
    // the band and tighter at its edges.
    const std::uint32_t N = 4 * f->taps, bins = N / 2;
    const auto betaFor = [](double db) {
        const double cap = std::pow(10.0, db / 20.0);
        return 1.0 / (4.0 * cap * cap);
    };
    const double betaMiddle = betaFor(p.max_boost_db), betaEdge = betaFor(p.edge_boost_db);
    std::vector<Matrix> C(bins + 1);
    std::vector<double> nearPower(bins + 1), weight(bins + 1);
    for (std::uint32_t k = 0; k <= bins; ++k) {
        const double hz = static_cast<double>(k) * rate / N;
        weight[k] = bandWeight(hz, p.low_hz, p.high_hz);
        const double s = edgeShape(hz, p.low_hz, p.high_hz);
        const double beta = std::pow(betaMiddle, 1.0 - s) * std::pow(betaEdge, s);
        const Complex reference = rawPlant(p, {0.0, 0.0, 0.0}, hz, c)[0][0];
        std::vector<Matrix> plants;
        Matrix normal{}, adjointSum{};
        for (const Head& head : heads) {
            Matrix h = rawPlant(p, head, hz, c);
            for (int i = 0; i < 2; ++i)
                for (int j = 0; j < 2; ++j) h[i][j] = weight[k] * h[i][j] / reference + (i == j ? 1.0 - weight[k] : 0.0);
            const Matrix ha = adjoint(h), hh = multiply(ha, h);
            for (int i = 0; i < 2; ++i)
                for (int j = 0; j < 2; ++j) {
                    normal[i][j] += hh[i][j];
                    adjointSum[i][j] += ha[i][j];
                }
            plants.push_back(h);
        }
        normal[0][0] += K * beta;
        normal[1][1] += K * beta;
        C[k] = multiply(inverse(normal), adjointSum);
        // What reaches each near ear from its own channel, over the heads: the canceller's own colour.
        double power = 0.0;
        for (const Matrix& h : plants) {
            const Matrix total = multiply(h, C[k]);
            power += 0.5 * (std::norm(total[0][0]) + std::norm(total[1][1]));
        }
        nearPower[k] = power / K;
    }

    // The correction: the inverse of that colour, smoothed over a third of an octave, held to 6 dB and only inside the band.
    if (p.equalize) {
        std::vector<double> prefix(bins + 2, 0.0);
        for (std::uint32_t k = 0; k <= bins; ++k) prefix[k + 1] = prefix[k] + nearPower[k];
        const double limit = std::pow(10.0, 6.0 / 20.0), third = std::pow(2.0, 1.0 / 6.0);
        for (std::uint32_t k = 1; k <= bins; ++k) {
            const auto lo = static_cast<std::uint32_t>(std::max(1.0, std::floor(k / third)));
            const auto hi = static_cast<std::uint32_t>(std::min<double>(bins, std::ceil(k * third)));
            const double mean = (prefix[hi + 1] - prefix[lo]) / (hi + 1 - lo);
            const double gain = std::clamp(1.0 / std::sqrt(std::max(mean, 1e-12)), 1.0 / limit, limit);
            const double applied = 1.0 + weight[k] * (gain - 1.0);
            for (auto& row : C[k])
                for (Complex& v : row) v *= applied;
        }
    }

    // Impulse responses by the inverse FFT, which sums without the 1/N; then a half-Hann fade over the last eighth of the taps.
    RealFft design(N);
    std::vector<float> re(bins), im(bins), time(N);
    for (int speaker = 0; speaker < 2; ++speaker)
        for (int channel = 0; channel < 2; ++channel) {
            const auto value = [&](std::uint32_t k) { return C[k][speaker][channel] * std::polar(1.0, -2.0 * kPiD * k / N * f->latencyFrames); };
            re[0] = static_cast<float>(value(0).real());
            im[0] = static_cast<float>(value(bins).real());
            for (std::uint32_t k = 1; k < bins; ++k) {
                const Complex v = value(k);
                re[k] = static_cast<float>(v.real());
                im[k] = static_cast<float>(v.imag());
            }
            design.inverse(re.data(), im.data(), time.data());
            std::vector<float>& out = f->impulse[speaker][channel];
            out.assign(f->taps, 0.0f);
            const std::uint32_t fade = f->taps / 8;
            for (std::uint32_t n = 0; n < f->taps; ++n) {
                double g = 1.0;
                if (n >= f->taps - fade) g = 0.5 * (1.0 + std::cos(kPiD * (n - (f->taps - fade) + 1) / fade));
                out[n] = static_cast<float>(time[n] / N * g);
            }
        }

    // Partition spectra, as the engine's convolutions take them: each block of taps zero-padded to two blocks.
    RealFft fft(2 * B);
    f->re.reset(4ull * f->partitions * B);
    f->im.reset(4ull * f->partitions * B);
    std::vector<float> segment(2 * B);
    for (std::uint32_t pair = 0; pair < 4; ++pair) {
        const std::vector<float>& taps = f->impulse[pair / 2][pair % 2];
        for (std::uint32_t part = 0; part < f->partitions; ++part) {
            std::fill(segment.begin(), segment.end(), 0.0f);
            for (std::uint32_t t = 0; t < B && part * B + t < f->taps; ++t) segment[t] = taps[part * B + t];
            const std::size_t at = (static_cast<std::size_t>(pair) * f->partitions + part) * B;
            fft.forward(segment.data(), f->re.data() + at, f->im.data() + at);
        }
    }
    return f;
}

}  // namespace mspatial::core
