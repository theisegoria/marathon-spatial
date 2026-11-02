// SPDX-License-Identifier: Apache-2.0
// The engine's two halves and the boundary between them.
//
// Control side -> render side: one snapshot of all continuous state, published by a wait-free triple index, and one ring of
// structural commands. Render side -> control side: a ring of retired objects, so nothing is freed on the render side, and
// relaxed atomics for telemetry. All render memory is allocated at creation.
#pragma once
#include <os/lock.h>

#include <memory>
#include <vector>

#include "Ambisonics.hpp"
#include "Common.hpp"
#include "Dsp.hpp"
#include "Fft.hpp"
#include "hrtf/HrtfSet.hpp"

namespace mspatial::core {

struct BufferObject {
    Aligned<float> data;  // Interleaved
    std::uint64_t frames{};
    std::uint32_t channels{};
};

struct SourceSnapshot {
    ms_source_params params;
    std::uint32_t reflectionCount;
    ms_reflection_tap reflections[MS_MAX_REFLECTION_TAPS];
};
struct BedSnapshot {
    Quat orientation;
    float gain;
    std::uint32_t headLocked;
    std::uint32_t bus;  // Index; 0 is the master.
};
struct BusSnapshot {
    float gain, highpass, lowpass;
};
struct ListenerSnapshot {
    Vec3 position, velocity, forward{0, 0, -1}, up{0, 1, 0};
};

enum class CommandType : std::uint32_t {
    SourceReset, SourceBuffer, SourceCallback, SourcePlay, SourceStop, SourceSeek,
    BedReset, BedBuffer, BedCallback, BedPlay, BedStop,
    BusProcessor, BlockCallback, BufferRetire, SetHrtf
};

struct Command {
    CommandType type;
    std::uint32_t index;
    std::uint32_t a;
    std::uint32_t b;
    std::uint64_t frame;
    double position;
    void* pointer;
    void* pointer2;
};

struct Retired {
    enum class Kind : std::uint32_t { Buffer, Hrtf } kind;
    void* pointer;
    void* pointer2;
};

// State of one convolved voice. Pooled: as many as max_hrtf_voices.
struct HrtfVoice {
    float* input;        // [2 B]
    float* fdlRe;        // [P][B]
    float* fdlIm;
    float* filterRe;     // [slot][ear][P][B]: the two most recent filters
    float* filterIm;
    float* delay[2];     // Each ear's line, written twice so a read window never wraps. [2 kDelayFrames]
    float* bands[3];     // Reflection delay lines, low, mid, high. [kReflectionFrames]
    std::uint32_t fdlHead;
    std::uint32_t current;      // Which of filterRe / filterIm holds the newest.
    std::uint32_t partitions[2];
    std::uint32_t delayWrite;
    std::uint32_t bandWrite;
    float delayNow[2];
    float delayTarget[2];
    Vec3 lastDirection[2];
    const hrtf::HrtfSet* lastSet;
    bool haveFilter;
    bool crossfade;
    Biquad near[2];
    float nearGain[2];
    Biquad splitLow, splitHigh;
    struct Tap {
        float delay;
        float gain[3];
        float encode[16];
        bool active;
    } taps[MS_MAX_REFLECTION_TAPS];
};

struct RenderSource {
    enum class State : std::uint32_t { Idle, Playing, Stopping, Tail } state;
    BufferObject* buffer;
    std::uint32_t channel;
    bool loop;
    ms_source_callback callback;
    void* user;
    double position;           // Frames into the content; negative before it starts.
    double seekTarget;
    std::uint64_t startFrame;  // Pending start while scheduled.
    bool scheduled;
    bool primed;
    Vec3 smoothed;
    Vec3 smoothedDirection;
    float gain;
    float dryGain;             // What the reflections hear: follows params.gain, and falls to zero on a stop.
    float rate;
    float pan[2];              // Last equal-power gains, MS_SPATIAL_PAN
    // Pull buffer for callback sources: older frames first, for the interpolator's history.
    float* pull;
    std::uint32_t pullFill;
    float pullFraction;
    bool pullEnded;
    bool contentEnded;
    std::uint32_t tailBlocks;
    enum Stage { HighPass, LowPass, EqLow, EqMid, EqHigh, Air1, Air2, Rear, StageCount };
    Biquad stages[StageCount];
    float stageValue[StageCount];
    HrtfVoice* voice;
    std::uint32_t voiceDrain;  // Blocks a demoted voice still rings out: filter tail, ear delay and reflection delay lines.
    float directMix;
    bool wantDirect;
    float score;
    float encode[kMaxChannels];
    float encodeWeight;
    float peak;
};

struct RenderBed {
    enum class State : std::uint32_t { Idle, Playing, Stopping } state;
    BufferObject* buffer;
    ms_bed_callback callback;
    void* user;
    std::uint32_t order;
    bool loop;
    bool scheduled;
    bool haveRotation;
    bool haveGain;
    std::uint64_t startFrame;
    std::int64_t position;
    float gain;
    float rotation[kRotationCells];
};

struct RenderBus {
    float* left;
    float* right;
    float gain;
    Biquad highpass[2], lowpass[2];
    float highpassHz, lowpassHz;
    ms_bus_processor processor;
    void* user;
    std::uint32_t parent;
    bool used;
    bool primed;  // A new bus starts at its gain instead of gliding there.
};

inline constexpr float kAirShelfHz[2] = {4000.0f, 12000.0f};
inline constexpr std::uint32_t kDelayFrames = 1024;  // Fits the largest onset delay, a block and the interpolator.
inline constexpr std::uint32_t kReflectionFrames = 8192;
inline constexpr std::uint32_t kMaxFilterTaps = 512;

struct AtomicStats {
    std::atomic<std::uint64_t> frames{0}, blocks{0}, overs{0}, faults{0};
    std::atomic<std::uint32_t> playing{0}, direct{0}, ambisonic{0};
    std::atomic<float> peak[2]{}, rms[2]{};
    std::atomic<double> last{0.0}, max{0.0};
    std::atomic<std::uint64_t> forward{0}, inverse{0}, spectrum{0}, real{0}, updates{0}, opBlocks{0};
};

struct SourcePublished {
    std::atomic<std::uint32_t> playing{0}, direct{0};
    std::atomic<double> position{0.0};
    std::atomic<float> level{0.0f};
};

}  // namespace mspatial::core

struct ms_engine {
    using Snapshotted = mspatial::core::TripleIndex;
    ms_engine_desc desc{};
    float rate{};
    double ticksToSeconds{};
    std::uint32_t B{}, C{}, maxPartitions{};
    float inverseScale{};
    float smoothing{};

    // Control side.
    os_unfair_lock lock = OS_UNFAIR_LOCK_INIT;
    struct Slot {
        bool alive{};
        std::uint16_t generation{1};
    };
    std::vector<Slot> sourceSlots, bedSlots, bufferSlots, busSlots;
    std::vector<std::unique_ptr<mspatial::core::BufferObject>> buffers, retiringBuffers;
    std::vector<mspatial::core::SourceSnapshot> workSources;
    std::vector<mspatial::core::BedSnapshot> workBeds;
    std::vector<mspatial::core::BusSnapshot> workBuses;
    mspatial::core::ListenerSnapshot workListener;
    std::vector<mspatial::core::Command> pending;
    std::size_t pendingHead{};
    mspatial::hrtf::HrtfSet* currentSet{};
    std::shared_ptr<const mspatial::core::AmbisonicDecoder> currentDecoder;
    struct Held {
        mspatial::hrtf::HrtfSet* set;
        std::shared_ptr<const mspatial::core::AmbisonicDecoder> decoder;
    };
    std::vector<Held> heldSets;
    std::uint32_t highestSource{};

    // The boundary.
    mspatial::core::TripleIndex triple;
    mspatial::core::Aligned<mspatial::core::SourceSnapshot> snapSources[3];
    mspatial::core::Aligned<mspatial::core::BedSnapshot> snapBeds[3];
    mspatial::core::Aligned<mspatial::core::BusSnapshot> snapBuses[3];
    mspatial::core::ListenerSnapshot snapListener[3];
    mspatial::core::SpscRing<mspatial::core::Command> commands;
    mspatial::core::SpscRing<mspatial::core::Retired> retired;
    std::atomic<std::uint64_t> headPacked{0};
    mspatial::core::AtomicStats stats;
    mspatial::core::Aligned<mspatial::core::SourcePublished> published;

    // Render side.
    mspatial::core::RealFft fft;
    mspatial::core::ShRotation rotation;
    std::uint64_t frame{};
    mspatial::core::Aligned<mspatial::core::RenderSource> sources;
    mspatial::core::Aligned<mspatial::core::RenderBed> beds;
    mspatial::core::Aligned<mspatial::core::RenderBus> buses;
    mspatial::core::Aligned<mspatial::core::HrtfVoice> voices;
    mspatial::core::Aligned<std::uint32_t> freeVoices;
    std::uint32_t freeVoiceCount{};
    mspatial::core::Aligned<float> voiceMemory, pullMemory, busMemory, scratch, shMix, shInput, shFdlRe, shFdlIm;
    mspatial::core::Aligned<float> crossfadeRamp;
    mspatial::core::Aligned<std::uint32_t> ranking;
    std::uint32_t shFdlHead{};
    std::uint32_t shActiveBlocks{};
    const mspatial::hrtf::HrtfSet* set{};
    const mspatial::core::AmbisonicDecoder* decoder{};
    const mspatial::core::AmbisonicDecoder* previousDecoder{};
    mspatial::hrtf::HrtfSet* retiringSet{};
    mspatial::core::Aligned<float> fifo;
    std::uint32_t fifoRead{}, fifoFill{};
    float limiterGain{1.0f};
    float airTable[2]{};   // Shelf gains in dB per metre
    std::uint64_t opForward{}, opInverse{}, opSpectrum{}, opReal{}, opUpdates{};
    mspatial::core::Quat headFrom;
    mspatial::core::Vec3 right{1, 0, 0}, up{0, 1, 0}, back{0, 0, 1};
    const mspatial::core::ListenerSnapshot* listener{};
    const mspatial::core::SourceSnapshot* frameSources{};
    const mspatial::core::BedSnapshot* frameBeds{};
    const mspatial::core::BusSnapshot* frameBuses{};
    // Gain from a bus's content to the master's input, leaving out the master's own gain: content that joins the master
    // directly (the ambisonic decode) goes through that gain once there. busProductFrom is the previous block's value.
    float busProduct[256]{};
    float busProductFrom[256]{};
    ms_block_callback blockCallback{};
    void* blockUser{};
};
