// SPDX-License-Identifier: Apache-2.0
// Real FFT on Accelerate in split form. Bin 0 holds DC in re and Nyquist in im, as vDSP packs them.
#pragma once
#include <Accelerate/Accelerate.h>

#include <cstdint>
#include <stdexcept>

namespace mspatial::core {

class RealFft {
public:
    RealFft() = default;
    explicit RealFft(std::uint32_t size) { init(size); }
    RealFft(const RealFft&) = delete;
    RealFft& operator=(const RealFft&) = delete;
    ~RealFft() {
        if (setup_) vDSP_destroy_fftsetup(setup_);
    }
    void init(std::uint32_t size) {
        if (setup_) vDSP_destroy_fftsetup(setup_);
        size_ = size;
        log2n_ = 0;
        while ((1u << log2n_) < size) ++log2n_;
        if ((1u << log2n_) != size) throw std::invalid_argument("FFT size");
        setup_ = vDSP_create_fftsetup(log2n_, kFFTRadix2);
        if (!setup_) throw std::bad_alloc();
    }
    std::uint32_t size() const noexcept { return size_; }
    std::uint32_t bins() const noexcept { return size_ / 2; }

    // re and im each hold size / 2 values. The result is twice the mathematical DFT.
    void forward(const float* time, float* re, float* im) const noexcept {
        DSPSplitComplex z{re, im};
        vDSP_ctoz(reinterpret_cast<const DSPComplex*>(time), 2, &z, 1, size_ / 2);
        vDSP_fft_zrip(setup_, &z, 1, log2n_, kFFTDirection_Forward);
    }
    // Destroys re and im. The result is size times the mathematical inverse, so forward then inverse gains 2 x size.
    void inverse(float* re, float* im, float* time) const noexcept {
        DSPSplitComplex z{re, im};
        vDSP_fft_zrip(setup_, &z, 1, log2n_, kFFTDirection_Inverse);
        vDSP_ztoc(&z, 1, reinterpret_cast<DSPComplex*>(time), 2, size_ / 2);
    }

private:
    FFTSetup setup_{};
    std::uint32_t size_{};
    std::uint32_t log2n_{};
};

}  // namespace mspatial::core
