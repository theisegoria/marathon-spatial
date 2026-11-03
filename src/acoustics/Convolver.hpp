// Uniformly partitioned overlap-save convolution of one input with a multichannel impulse response, on vDSP.
//
// The response is cut into blocks of `block` samples, each held as a spectrum; the input's spectra sit in a frequency-domain delay line,
// so every block costs one forward FFT plus, per output channel, one multiply-accumulate per partition and one inverse FFT. Latency is
// exactly one block. A new response is handed over without waiting: the control thread fills a free slot and publishes it, and the render
// thread picks it up at a block boundary and crossfades from the old response over that block, so a changing room never clicks.
//
// prepare() allocates. stage() never allocates. process() never allocates, locks or makes a system call.
#pragma once

#include "core/Fft.hpp"

#include <atomic>
#include <cstdint>
#include <vector>

namespace mspatial::acoustics {

class PartitionedConvolver {
public:
    // block must be a power of two. maxSamples bounds every response this object will take.
    void prepare(std::uint32_t block, std::uint32_t maxSamples, std::uint32_t channels);
    bool prepared() const { return block_ > 0; }

    std::uint32_t block() const { return block_; }
    std::uint32_t channels() const { return channels_; }
    std::uint32_t latency() const { return block_; }

    // Control thread. ir[c] holds `length` samples for channel c; length is cut to maxSamples. False when the previous response has
    // not been picked up yet, in which case nothing changes and the caller tries again later.
    bool stage(const float* const* ir, std::uint32_t length);
    // Control thread. A silent response: the output fades to zero over one block.
    bool stageSilence();

    // Render thread. outputs[c] receives frames samples for each channel; a null input means silence.
    void process(const float* input, float* const* outputs, std::uint32_t frames);
    // Render thread. Forgets the input history.
    void reset();

private:
    static constexpr std::uint32_t kSlots = 3;
    static constexpr std::uint32_t kNoSlot = 0xFFu;

    void runBlock();
    void accumulate(std::uint32_t slot, std::uint32_t channel, float* re, float* im) const;
    float* spectrum(std::uint32_t slot, std::uint32_t channel, std::uint32_t partition) {
        return slots_.data() + ((std::size_t(slot) * channels_ + channel) * partitions_ + partition) * std::size_t(2) * block_;
    }
    const float* spectrum(std::uint32_t slot, std::uint32_t channel, std::uint32_t partition) const {
        return slots_.data() + ((std::size_t(slot) * channels_ + channel) * partitions_ + partition) * std::size_t(2) * block_;
    }

    core::RealFft fft_;
    std::uint32_t block_ = 0;
    std::uint32_t channels_ = 0;
    std::uint32_t partitions_ = 0;
    std::vector<float> slots_;                       // [slot][channel][partition] spectra: block re then block im.
    std::uint32_t slotPartitions_[kSlots] = {};      // Partitions in use per slot.
    std::vector<float> history_;                     // Input spectra, a ring of `partitions_`, block re then block im each.
    std::uint32_t head_ = 0;
    std::vector<float> window_;                      // The last two input blocks, time domain.
    std::vector<float> scratch_;                     // FFT work: two blocks of re and im each.
    std::vector<float> staging_;                     // Control thread's time-domain partition.
    std::vector<float> inFifo_;
    std::vector<float> outFifo_;                     // [channel][block], the block being played.
    std::uint32_t fill_ = 0;                         // Samples of the current input block received; also the read position in outFifo_.

    std::atomic<std::uint32_t> pending_{kNoSlot};    // Written by stage, taken by the render thread.
    std::atomic<std::uint32_t> inUse_{0};            // Bit mask of slots the render thread may read.
    std::uint32_t active_ = kNoSlot;
    std::uint32_t retiring_ = kNoSlot;
};

}  // namespace mspatial::acoustics
