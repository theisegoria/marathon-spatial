// Late reverberation, one feedback delay network per room, exchanging energy through portals.
//
// Each room is a 16-line lossless Householder network. The per-line attenuation follows the room's octave-band decay rates exactly as the
// simulator measured them, so T60 per band comes from real volume and absorption. Levels are physical: the output is pressure referenced to the
// free-field direct sound at 1 m, so the late energy of an impulse follows 16 pi / A with no tuned wet level. Portals carry energy between
// networks at the diffuse-field rate c S / 4 V, which reproduces the double-slope decay of coupled rooms.
//
// process() is real time: no locks, no allocation, no system calls. apply() is too, and may run every audio block.
#pragma once

#include "acoustics/AcousticSnapshot.hpp"

#include <array>
#include <atomic>
#include <vector>

namespace mspatial::acoustics {

inline constexpr std::uint32_t kReverbLines = 16;
inline constexpr std::uint32_t kReverbSections = kBandCount;
inline constexpr std::uint32_t kReverbChunk = 32;
// Binaural output: six bands split at these frequencies (Linkwitz-Riley), each with its own interaural coherence. Band centres are the
// geometric means, the outer two an octave beyond: 141, 283, 566, 1131, 2263 and 4525 Hz.
inline constexpr std::uint32_t kCoherenceBands = 6;
inline constexpr std::array<float, kCoherenceBands - 1> kCoherenceCrossovers{200.f, 400.f, 800.f, 1600.f, 3200.f};

struct LateReverbConfig {
    float sampleRate = 48000.f;
    float speedOfSound = kDefaultSpeedOfSound;
    // Mean line delay per room in seconds; empty entries or absent rooms use 30 ms. Fixed for the life of prepare(): a room whose delays changed
    // live would change pitch. A good value is meanFreePath / c.
    std::vector<float> meanDelaySeconds;
    std::uint32_t rooms = kMaxRooms;
    std::uint32_t links = kMaxRoomLinks;
    // Longest delay setHybrid may ask for, in samples. 0 means no hybrid room.
    std::uint32_t maxHybridDelay = 0;
};

class LateReverbNetwork {
public:
    // Allocates. Not real time.
    void prepare(const LateReverbConfig& config);

    // Takes new room and link parameters. Decay changes glide over a few milliseconds.
    void apply(const AcousticSnapshot& snapshot);

    // inputs[r] is the mono send into room r, or null. outputs hold left and right for room r, or null. Rooms with no active state output silence.
    void process(const float* const* inputs, float* const* outputsLeft, float* const* outputsRight, std::uint32_t frames);

    void reset();

    // Hybrid with a convolution reverb: room's own input reaches its network `delaySamples` later and scaled per band by gainDb, so the
    // convolution owns the start of the tail and the network carries on from where the response fades. Inputs injected through portals
    // are not delayed. kNone turns it off. Real time: call from the render thread before process(); a changed room clears the delay.
    void setHybrid(std::uint32_t room, std::uint32_t delaySamples, const Bands& gainDb);

    // Binaural output: left and right become a.A + b.B and a.A - b.B of the network's two orthogonal outputs, with a^2 + b^2 = 1 and
    // (a^2 - b^2) the given interaural coherence per band, so the tail reaches the ears as a diffuse field would through the listener's
    // HRTF set, with its level unchanged. Null returns to the two orthogonal outputs. Any thread; the render picks it up at its next chunk.
    void setCoherence(const float* coherence);

    // Delay memory in one room's network, in samples.
    std::uint32_t totalDelay(std::uint32_t room) const;

private:
    struct Biquad {
        float b0 = 1.f, b1 = 0.f, b2 = 0.f, a1 = 0.f, a2 = 0.f;
    };
    struct Filter {
        std::array<Biquad, kReverbSections> sections;
        std::array<float, kReverbSections> z1{}, z2{};
    };
    struct Room {
        bool active = false;
        bool designed = false;
        std::array<std::uint32_t, kReverbLines> lengths{};
        std::array<std::uint32_t, kReverbLines> positions{};
        std::array<std::uint32_t, kReverbLines> offsets{};
        std::array<float, kReverbLines> flatGain{};
        std::array<Filter, kReverbLines> filters;
        Bands targetDbPerSample{};
        Bands currentDbPerSample{};
        float totalDelay = 0.f;
        float inputGain = 0.f;
        std::uint32_t silentChunks = 0;
    };
    struct Link {
        bool active = false;
        std::uint32_t from = kNone, to = kNone;
        float gain = 0.f;
        bool shaped = false;
        float phase = 0.f;
        float phaseStep = 0.f;
        Filter filter;
    };

    void design(Room& room);
    void designLink(Link& link, const Bands& rate, float fromVolume, float toVolume);
    void renderRoom(std::uint32_t index, const float* input, std::uint32_t offset, std::uint32_t frames, float* left, float* right);
    // The binaural band split of one signal: per crossover, a low-pass and a high-pass of two Butterworth stages each.
    struct Section {float z1 = 0.f, z2 = 0.f;};
    struct Split {std::array<std::array<Section, 4>, kCoherenceBands - 1> stages{};};
    void binaural(std::uint32_t room, float* left, float* right, std::uint32_t frames);
    std::array<Biquad, kCoherenceBands - 1> splitLow_{}, splitHigh_{};
    std::vector<std::array<Split, 2>> splits_;
    std::array<std::atomic<float>, kCoherenceBands> coherence_{};
    std::atomic<bool> binaural_{false};

    LateReverbConfig config_;
    std::vector<Room> rooms_;
    std::vector<Link> links_;
    std::vector<float> delayMemory_;
    std::vector<std::array<float, kReverbChunk>> roomSend_, roomSendOther_;
    std::vector<std::array<float, kReverbChunk>> linkInput_;
    std::array<std::array<float, kBandCount>, kBandCount> inverseResponse_{};

    struct Hybrid {
        std::uint32_t room = kNone;
        std::uint32_t delay = 0;
        std::uint32_t write = 0;
        Bands gainDb{};
        float flatGain = 1.f;
        bool shaped = false;
        Filter filter;
        std::array<float, kReverbChunk> chunk{};
    } hybrid_;
    std::vector<float> hybridMemory_;
};

}  // namespace mspatial::acoustics
