#include "acoustics/Convolver.hpp"

#include <algorithm>
#include <cstring>
#include <stdexcept>

namespace mspatial::acoustics {

void PartitionedConvolver::prepare(std::uint32_t block, std::uint32_t maxSamples, std::uint32_t channels) {
    if (block < 16 || (block & (block - 1)) != 0 || channels == 0) throw std::invalid_argument("convolver block");
    fft_.init(block * 2);
    block_ = block;
    channels_ = channels;
    partitions_ = std::max<std::uint32_t>(1, (maxSamples + block - 1) / block);
    slots_.assign(std::size_t(kSlots) * channels * partitions_ * 2 * block, 0.f);
    std::fill(std::begin(slotPartitions_), std::end(slotPartitions_), 0u);
    history_.assign(std::size_t(partitions_) * 2 * block, 0.f);
    window_.assign(std::size_t(2) * block, 0.f);
    scratch_.assign(std::size_t(5) * block, 0.f);
    staging_.assign(std::size_t(2) * block, 0.f);
    inFifo_.assign(block, 0.f);
    outFifo_.assign(std::size_t(channels) * block, 0.f);
    head_ = 0;
    fill_ = 0;
    pending_.store(kNoSlot, std::memory_order_relaxed);
    inUse_.store(0, std::memory_order_relaxed);
    active_ = kNoSlot;
    retiring_ = kNoSlot;
}

bool PartitionedConvolver::stage(const float* const* ir, std::uint32_t length) {
    if (!prepared() || pending_.load(std::memory_order_acquire) != kNoSlot) return false;
    const std::uint32_t busy = inUse_.load(std::memory_order_acquire);
    std::uint32_t slot = 0;
    while (slot < kSlots && (busy & (1u << slot))) ++slot;
    if (slot == kSlots) return false;
    length = std::min(length, partitions_ * block_);
    const std::uint32_t used = (length + block_ - 1) / block_;
    // forward() and inverse() together gain 2N per transform pair, and both input and response pass forward(): 4N in all.
    const float scale = 1.f / (4.f * float(2 * block_));
    for (std::uint32_t c = 0; c < channels_; ++c) {
        for (std::uint32_t k = 0; k < used; ++k) {
            std::fill(staging_.begin(), staging_.end(), 0.f);
            const std::uint32_t first = k * block_, count = std::min(block_, length - first);
            if (ir && ir[c]) std::memcpy(staging_.data(), ir[c] + first, count * sizeof(float));
            float* re = spectrum(slot, c, k);
            float* im = re + block_;
            fft_.forward(staging_.data(), re, im);
            vDSP_vsmul(re, 1, &scale, re, 1, std::size_t(2) * block_);
        }
    }
    slotPartitions_[slot] = used;
    pending_.store(slot, std::memory_order_release);
    return true;
}

bool PartitionedConvolver::stageSilence() { return stage(nullptr, 0); }

void PartitionedConvolver::reset() {
    std::fill(history_.begin(), history_.end(), 0.f);
    std::fill(window_.begin(), window_.end(), 0.f);
    std::fill(inFifo_.begin(), inFifo_.end(), 0.f);
    std::fill(outFifo_.begin(), outFifo_.end(), 0.f);
    fill_ = 0;
}

void PartitionedConvolver::accumulate(std::uint32_t slot, std::uint32_t channel, float* re, float* im) const {
    std::fill_n(re, block_, 0.f);
    std::fill_n(im, block_, 0.f);
    if (slot == kNoSlot) return;
    const std::uint32_t used = slotPartitions_[slot];
    DSPSplitComplex acc{re + 1, im + 1};
    for (std::uint32_t k = 0; k < used; ++k) {
        const float* x = history_.data() + std::size_t((head_ + partitions_ - k) % partitions_) * 2 * block_;
        const float* h = spectrum(slot, channel, k);
        // Bin 0 packs DC in re and Nyquist in im, both real.
        re[0] += x[0] * h[0];
        im[0] += x[block_] * h[block_];
        DSPSplitComplex xs{const_cast<float*>(x) + 1, const_cast<float*>(x) + block_ + 1};
        DSPSplitComplex hs{const_cast<float*>(h) + 1, const_cast<float*>(h) + block_ + 1};
        vDSP_zvma(&xs, 1, &hs, 1, &acc, 1, &acc, 1, block_ - 1);
    }
}

void PartitionedConvolver::runBlock() {
    bool crossfade = false;
    const std::uint32_t pending = pending_.load(std::memory_order_acquire);
    if (pending != kNoSlot) {
        retiring_ = active_;
        active_ = pending;
        inUse_.store((1u << active_) | (retiring_ != kNoSlot ? 1u << retiring_ : 0u), std::memory_order_release);
        pending_.store(kNoSlot, std::memory_order_release);
        crossfade = true;
    }

    // Slide the input window and take the newest spectrum into the delay line.
    std::memmove(window_.data(), window_.data() + block_, block_ * sizeof(float));
    std::memcpy(window_.data() + block_, inFifo_.data(), block_ * sizeof(float));
    head_ = (head_ + 1) % partitions_;
    float* x = history_.data() + std::size_t(head_) * 2 * block_;
    fft_.forward(window_.data(), x, x + block_);

    float* re = scratch_.data();
    float* im = re + block_;
    float* time = im + block_;          // Two blocks.
    float* old = time + 2 * block_;     // One block.
    for (std::uint32_t c = 0; c < channels_; ++c) {
        float* out = outFifo_.data() + std::size_t(c) * block_;
        if (crossfade) {
            accumulate(retiring_, c, re, im);
            fft_.inverse(re, im, time);
            std::memcpy(old, time + block_, block_ * sizeof(float));
        }
        accumulate(active_, c, re, im);
        fft_.inverse(re, im, time);
        // Overlap-save keeps the second half: the first is wrapped by the circular convolution.
        std::memcpy(out, time + block_, block_ * sizeof(float));
        if (crossfade) {
            const float step = 1.f / float(block_);
            for (std::uint32_t n = 0; n < block_; ++n) {
                const float w = (float(n) + 0.5f) * step;
                out[n] = old[n] + (out[n] - old[n]) * w;
            }
        }
    }
    if (crossfade && retiring_ != kNoSlot) {
        retiring_ = kNoSlot;
        inUse_.store(1u << active_, std::memory_order_release);
    }
}

void PartitionedConvolver::process(const float* input, float* const* outputs, std::uint32_t frames) {
    std::uint32_t done = 0;
    while (done < frames) {
        const std::uint32_t span = std::min(frames - done, block_ - fill_);
        if (input)
            std::memcpy(inFifo_.data() + fill_, input + done, span * sizeof(float));
        else
            std::fill_n(inFifo_.data() + fill_, span, 0.f);
        for (std::uint32_t c = 0; c < channels_; ++c)
            if (outputs && outputs[c]) std::memcpy(outputs[c] + done, outFifo_.data() + std::size_t(c) * block_ + fill_, span * sizeof(float));
        fill_ += span;
        done += span;
        if (fill_ == block_) {
            runBlock();
            fill_ = 0;
        }
    }
}

}  // namespace mspatial::acoustics
