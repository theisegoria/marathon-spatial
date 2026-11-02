// SPDX-License-Identifier: Apache-2.0
#define ACCELERATE_NEW_LAPACK
#include "Ambisonics.hpp"

#include <dispatch/dispatch.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <vector>

#include "Fft.hpp"
#include "hrtf/HrtfSet.hpp"

namespace mspatial::core {

namespace {

struct Sn3dTable {
    double norm[kMaxOrder + 1][kMaxOrder + 1];
    Sn3dTable() {
        for (std::uint32_t l = 0; l <= kMaxOrder; ++l)
            for (std::uint32_t m = 0; m <= l; ++m) {
                double ratio = 1.0;
                for (std::uint32_t k = l - m + 1; k <= l + m; ++k) ratio /= static_cast<double>(k);
                norm[l][m] = std::sqrt((m == 0 ? 1.0 : 2.0) * ratio);
            }
    }
};
const Sn3dTable gSn3d;

}  // namespace

void encodeDirection(Vec3 d, std::uint32_t order, float* gains) noexcept {
    const double x = d.x, y = d.y, z = d.z;
    double c = 1.0, s = 0.0;       // cos^m(el) cos(m az), cos^m(el) sin(m az)
    double pmm = 1.0;              // Double factorial (2m - 1)!!
    for (std::uint32_t m = 0; m <= order; ++m) {
        double p0 = pmm;
        double p1 = 0.0;
        for (std::uint32_t l = m; l <= order; ++l) {
            const double value = gSn3d.norm[l][m] * p0;
            if (m == 0) {
                gains[l * l + l] = static_cast<float>(value);
            } else {
                gains[l * l + l + m] = static_cast<float>(value * c);
                gains[l * l + l - m] = static_cast<float>(value * s);
            }
            const double next = l == m ? z * (2.0 * m + 1.0) * p0
                                       : (z * (2.0 * l + 1.0) * p0 - (l + m) * p1) / static_cast<double>(l + 1 - m);
            p1 = p0;
            p0 = next;
        }
        const double cc = c * x - s * y;
        s = s * x + c * y;
        c = cc;
        pmm *= 2.0 * m + 1.0;
    }
}

ShRotation::ShRotation() {
    const int L = static_cast<int>(kMaxOrder);
    table_.reset((kMaxOrder + 1ull) * (2 * kMaxOrder + 1) * (2 * kMaxOrder + 1));
    for (int l = 1; l <= L; ++l)
        for (int a = -l; a <= l; ++a)
            for (int b = -l; b <= l; ++b) {
                const double d = a == 0 ? 1.0 : 0.0;
                const double denominator = std::abs(b) == l ? 2.0 * l * (2.0 * l - 1.0) : static_cast<double>(l * l - b * b);
                const int aa = std::abs(a);
                table_[index(l, a, b)] = {std::sqrt((l * l - a * a) / denominator),
                                          0.5 * std::sqrt((1.0 + d) * (l + aa - 1) * (l + aa) / denominator) * (1.0 - 2.0 * d),
                                          -0.5 * std::sqrt(std::max(0, (l - aa - 1) * (l - aa)) / denominator) * (1.0 - d)};
            }
}

void ShRotation::compute(const float m3[9], std::uint32_t order, float* out) const noexcept {
    // Degree 1 is the rotation itself in the channel order Y, Z, X.
    static constexpr int kAxis[3] = {1, 2, 0};
    double r1[3][3];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) r1[i][j] = m3[kAxis[i] * 3 + kAxis[j]];
    const auto R = [&](int i, int j) { return r1[i + 1][j + 1]; };
    double blocks[kMaxOrder + 1][(2 * kMaxOrder + 1) * (2 * kMaxOrder + 1)];
    const auto M = [&](int l, int a, int b) -> double& { return blocks[l][(a + l) * (2 * l + 1) + (b + l)]; };
    M(0, 0, 0) = 1.0;
    for (int a = -1; a <= 1; ++a)
        for (int b = -1; b <= 1; ++b) M(1, a, b) = R(a, b);
    for (int l = 2; l <= static_cast<int>(order); ++l) {
        const auto P = [&](int i, int a, int b) {
            if (b == l) return R(i, 1) * M(l - 1, a, l - 1) - R(i, -1) * M(l - 1, a, -l + 1);
            if (b == -l) return R(i, 1) * M(l - 1, a, -l + 1) + R(i, -1) * M(l - 1, a, l - 1);
            return R(i, 0) * M(l - 1, a, b);
        };
        for (int a = -l; a <= l; ++a)
            for (int b = -l; b <= l; ++b) {
                const Coefficients& c = table_[index(l, a, b)];
                double sum = 0.0;
                if (c.u != 0.0) sum += c.u * P(0, a, b);
                if (c.v != 0.0) {
                    double v;
                    if (a == 0) v = P(1, 1, b) + P(-1, -1, b);
                    else if (a > 0) v = P(1, a - 1, b) * (a == 1 ? M_SQRT2 : 1.0) - P(-1, -a + 1, b) * (a == 1 ? 0.0 : 1.0);
                    else v = P(1, a + 1, b) * (a == -1 ? 0.0 : 1.0) + P(-1, -a - 1, b) * (a == -1 ? M_SQRT2 : 1.0);
                    sum += c.v * v;
                }
                if (c.w != 0.0) {
                    const double w = a > 0 ? P(1, a + 1, b) + P(-1, -a - 1, b) : P(1, a - 1, b) - P(-1, -a + 1, b);
                    sum += c.w * w;
                }
                M(l, a, b) = sum;
            }
    }
    std::size_t cursor = 0;
    for (int l = 0; l <= static_cast<int>(order); ++l)
        for (int a = -l; a <= l; ++a)
            for (int b = -l; b <= l; ++b) out[cursor++] = static_cast<float>(M(l, a, b));
}

std::unique_ptr<AmbisonicDecoder> buildDecoder(const hrtf::HrtfSet& set, std::uint32_t order, std::uint32_t blockFrames) {
    const std::uint32_t C = (order + 1) * (order + 1);
    const std::uint32_t Q = set.directions;
    if (Q < C) return nullptr;
    const int cInt = static_cast<int>(C), qInt = static_cast<int>(Q);

    // Enough length for filters plus spread delays, centred so a delay earlier than the mean stays causal.
    const auto spread = static_cast<std::uint32_t>(std::ceil(set.maxDelay));
    std::uint32_t N = 256;
    while (N < set.taps + 2 * spread + 32) N *= 2;
    const std::uint32_t bins = N / 2 + 1;

    double meanDelay = 0.0;
    for (std::size_t i = 0; i < Q * 2ull; ++i) meanDelay += set.delays[i];
    meanDelay /= Q * 2.0;
    const double shift = spread * 0.5 + 16.0;

    // The HRTFs as complex spectra with their delays restored: [ear][bin][direction].
    std::vector<float> hRe(static_cast<std::size_t>(Q) * bins * 2), hIm(hRe.size());
    float* hReal = hRe.data();
    float* hImag = hIm.data();
    const hrtf::HrtfSet* owner = &set;
    dispatch_apply(Q, DISPATCH_APPLY_AUTO, ^(std::size_t d) {
        RealFft fft(N);
        std::vector<float> time(N), re(N / 2), im(N / 2);
        for (std::uint32_t ear = 0; ear < 2; ++ear) {
            std::fill(time.begin(), time.end(), 0.0f);
            const float* filter = owner->filters.data() + (d * 2 + ear) * owner->taps;
            std::copy(filter, filter + owner->taps, time.begin());
            fft.forward(time.data(), re.data(), im.data());
            const double tau = owner->delays[d * 2 + ear] - meanDelay + shift;
            for (std::uint32_t k = 0; k < bins; ++k) {
                const std::complex<double> x = k == 0 ? std::complex<double>(re[0] * 0.5, 0.0)
                                                      : k == N / 2 ? std::complex<double>(im[0] * 0.5, 0.0)
                                                                   : std::complex<double>(re[k] * 0.5, im[k] * 0.5);
                const std::complex<double> y = x * std::polar(1.0, -2.0 * M_PI * k * tau / N);
                const std::size_t at = ((ear * static_cast<std::size_t>(bins)) + k) * Q + d;
                hReal[at] = static_cast<float>(y.real());
                hImag[at] = static_cast<float>(y.imag());
            }
        }
    });

    // Y: direction by channel. Solve (Y^T W Y + lambda I) X = Y^T W, then X maps a field sampled at the directions to channels.
    std::vector<float> Y(static_cast<std::size_t>(Q) * C);  // Row major.
    std::vector<float> weighted(Y.size());
    for (std::uint32_t q = 0; q < Q; ++q) {
        const Vec3 head{set.dirs[q * 3ull], set.dirs[q * 3ull + 1], set.dirs[q * 3ull + 2]};
        encodeDirection(headToAmbisonic(head), order, Y.data() + static_cast<std::size_t>(q) * C);
        const float weight = set.quadrature[q] * static_cast<float>(Q);
        for (std::uint32_t c = 0; c < C; ++c) weighted[static_cast<std::size_t>(q) * C + c] = Y[static_cast<std::size_t>(q) * C + c] * weight;
    }
    std::vector<double> gram(static_cast<std::size_t>(C) * C, 0.0);
    {
        std::vector<float> floatGram(gram.size());
        cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, cInt, cInt, qInt, 1.0f, weighted.data(), cInt, Y.data(), cInt, 0.0f, floatGram.data(), cInt);
        for (std::size_t i = 0; i < gram.size(); ++i) gram[i] = floatGram[i];
    }
    double trace = 0.0;
    for (std::uint32_t c = 0; c < C; ++c) trace += gram[c * (C + 1ull)];
    for (std::uint32_t c = 0; c < C; ++c) gram[c * (C + 1ull)] += 1e-4 * trace / C;
    std::vector<double> rhs(static_cast<std::size_t>(C) * Q);
    for (std::uint32_t q = 0; q < Q; ++q)
        for (std::uint32_t c = 0; c < C; ++c) rhs[static_cast<std::size_t>(q) * C + c] = weighted[static_cast<std::size_t>(q) * C + c];
    {
        char uplo = 'U';
        __LAPACK_int n = cInt, nrhs = qInt, lda = cInt, ldb = cInt, info = 0;
        dposv_(&uplo, &n, &nrhs, gram.data(), &lda, rhs.data(), &ldb, &info);
        if (info != 0) return nullptr;
    }
    std::vector<float> pinv(rhs.size());  // [direction][channel], i.e. X transposed: X(c, q) = pinv[q * C + c]
    for (std::size_t i = 0; i < rhs.size(); ++i) pinv[i] = static_cast<float>(rhs[i]);

    // cRe[ear][bin][channel]
    std::vector<float> cRe(static_cast<std::size_t>(bins) * C * 2), cIm(cRe.size());
    const double cutoffHz = static_cast<double>(order) * 343.0 / (2.0 * M_PI * set.headRadius);
    const auto cutoffBin = static_cast<std::uint32_t>(cutoffHz / (set.sampleRate / N));
    for (std::uint32_t ear = 0; ear < 2; ++ear) {
        std::vector<float> synthRe(Q), synthIm(Q), targetRe(Q), targetIm(Q);
        for (std::uint32_t k = 0; k < bins; ++k) {
            const float* re = hRe.data() + (static_cast<std::size_t>(ear) * bins + k) * Q;
            const float* im = hIm.data() + (static_cast<std::size_t>(ear) * bins + k) * Q;
            float* outRe = cRe.data() + (static_cast<std::size_t>(ear) * bins + k) * C;
            float* outIm = cIm.data() + (static_cast<std::size_t>(ear) * bins + k) * C;
            if (k > cutoffBin && k > 0) {
                // The phase the previous bin's solution produces, carried over and advanced by the common delay so the band
                // keeps its place in time; the measured magnitude.
                const float* lastRe = outRe - C;
                const float* lastIm = outIm - C;
                cblas_sgemv(CblasRowMajor, CblasNoTrans, qInt, cInt, 1.0f, Y.data(), cInt, lastRe, 1, 0.0f, synthRe.data(), 1);
                cblas_sgemv(CblasRowMajor, CblasNoTrans, qInt, cInt, 1.0f, Y.data(), cInt, lastIm, 1, 0.0f, synthIm.data(), 1);
                const float stepRe = static_cast<float>(std::cos(2.0 * M_PI * shift / N)), stepIm = static_cast<float>(-std::sin(2.0 * M_PI * shift / N));
                for (std::uint32_t q = 0; q < Q; ++q) {
                    const float magnitude = std::hypot(re[q], im[q]);
                    const float scale = magnitude / (std::hypot(synthRe[q], synthIm[q]) + 1e-12f);
                    targetRe[q] = scale * (synthRe[q] * stepRe - synthIm[q] * stepIm);
                    targetIm[q] = scale * (synthRe[q] * stepIm + synthIm[q] * stepRe);
                }
                cblas_sgemv(CblasRowMajor, CblasTrans, qInt, cInt, 1.0f, pinv.data(), cInt, targetRe.data(), 1, 0.0f, outRe, 1);
                cblas_sgemv(CblasRowMajor, CblasTrans, qInt, cInt, 1.0f, pinv.data(), cInt, targetIm.data(), 1, 0.0f, outIm, 1);
            } else {
                cblas_sgemv(CblasRowMajor, CblasTrans, qInt, cInt, 1.0f, pinv.data(), cInt, re, 1, 0.0f, outRe, 1);
                cblas_sgemv(CblasRowMajor, CblasTrans, qInt, cInt, 1.0f, pinv.data(), cInt, im, 1, 0.0f, outIm, 1);
            }
        }
    }

    auto decoder = std::make_unique<AmbisonicDecoder>();
    decoder->order = order;
    decoder->channels = C;
    decoder->blockFrames = blockFrames;
    decoder->partitions = (N + blockFrames - 1) / blockFrames;
    decoder->latencyFrames = static_cast<std::uint32_t>(std::lround(shift - meanDelay));
    const std::size_t total = 2ull * C * decoder->partitions * blockFrames;
    decoder->re.reset(total);
    decoder->im.reset(total);
    RealFft inverse(N), forward(blockFrames * 2);
    std::vector<float> time(N), re(N / 2), im(N / 2), padded(blockFrames * 2ull);
    for (std::uint32_t ear = 0; ear < 2; ++ear)
        for (std::uint32_t c = 0; c < C; ++c) {
            const float* inRe = cRe.data() + c + static_cast<std::size_t>(ear) * bins * C;
            const float* inIm = cIm.data() + c + static_cast<std::size_t>(ear) * bins * C;
            re[0] = inRe[0];
            im[0] = inRe[(N / 2) * C];
            for (std::uint32_t k = 1; k < N / 2; ++k) {
                re[k] = inRe[k * C];
                im[k] = inIm[k * C];
            }
            inverse.inverse(re.data(), im.data(), time.data());
            const float scale = 1.0f / static_cast<float>(N);
            const std::uint32_t fade = N / 16;
            for (std::uint32_t n = 0; n < N; ++n) {
                time[n] *= scale;
                if (n + fade >= N) {
                    const float t = static_cast<float>(N - n) / static_cast<float>(fade);
                    time[n] *= std::sin(0.5f * kPi * t) * std::sin(0.5f * kPi * t);
                }
            }
            for (std::uint32_t p = 0; p < decoder->partitions; ++p) {
                std::fill(padded.begin(), padded.end(), 0.0f);
                const std::uint32_t begin = p * blockFrames;
                for (std::uint32_t n = 0; n < blockFrames && begin + n < N; ++n) padded[n] = time[begin + n];
                forward.forward(padded.data(), const_cast<float*>(decoder->spectrumRe(ear, c, p)), const_cast<float*>(decoder->spectrumIm(ear, c, p)));
            }
        }
    return decoder;
}

}  // namespace mspatial::core

namespace mspatial::hrtf {

std::shared_ptr<const core::AmbisonicDecoder> HrtfSet::decoder(std::uint32_t order, std::uint32_t blockFrames) {
    const std::lock_guard<std::mutex> lock(decoderMutex_);
    for (const CachedDecoder& cached : decoders_)
        if (cached.order == order && cached.blockFrames == blockFrames) return cached.decoder;
    std::shared_ptr<const core::AmbisonicDecoder> built = core::buildDecoder(*this, order, blockFrames);
    if (built) decoders_.push_back({order, blockFrames, built});
    return built;
}

}  // namespace mspatial::hrtf
