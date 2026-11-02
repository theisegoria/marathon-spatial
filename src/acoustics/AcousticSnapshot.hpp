// What the audio thread reads: fixed-size plain data, published through a wait-free triple buffer. The audio thread never waits for
// the simulator or the GPU and never allocates; if nothing new arrived it keeps the last snapshot.
#pragma once

#include "acoustics/AcousticTypes.hpp"

#include <atomic>
#include <memory>

namespace mspatial::acoustics {

inline constexpr std::uint32_t kMaxSources = 64;
inline constexpr std::uint32_t kMaxPathsPerSource = 24;
inline constexpr std::uint32_t kMaxPortalPaths = 3;
inline constexpr std::uint32_t kMaxRooms = 32;
inline constexpr std::uint32_t kMaxRoomLinks = 64;

// One specular early reflection. A slot keeps its path for the path's whole life, so a renderer can bind one delay line to a slot and
// crossfade on id changes. gain is amplitude per band relative to the free-field direct sound at 1 m.
struct ReflectionPath {
    std::uint64_t id = 0;  // 0 for an empty slot. Stable across ticks for the same sequence of reflecting surfaces.
    float delaySeconds = 0.f;   // Whole path, source to listener: length / c.
    float excessSeconds = 0.f;  // After the direct sound: (length - direct distance) / c. This is what a tap delay is.
    Bands gain{};
    Vec3 arrival = vec3(0, 0, 0);   // Unit vector from the listener toward the last reflection.
    Vec3 emission = vec3(0, 0, 0);  // Unit vector from the source toward the first reflection.
    std::uint32_t order = 0;
};

// Sound through a chain of portals. arrival points at the portal nearest the listener, which is where the sound seems to come from.
struct PortalPath {
    Bands gain{};
    float delaySeconds = 0.f;   // Whole route / c.
    float excessSeconds = 0.f;  // After the straight line: (route - direct distance) / c.
    Vec3 arrival = vec3(0, 0, 0);
    Vec3 emission = vec3(0, 0, 0);
    std::uint32_t portals = 0;  // 0 for an empty entry.
};

struct SourceAcoustics {
    bool active = false;
    std::uint32_t room = kNone;
    float distance = 0.f;
    float visibility = 1.f;     // Fraction of the source disc with a clear line to the listener.
    Bands directGain{};         // Amplitude through whatever stands in the way, 1 when clear. Distance law not included.
    ReflectionPath reflections[kMaxPathsPerSource];
    PortalPath portalPaths[kMaxPortalPaths];
};

// Late-field state of one room. decayRate is the energy decay rate per second, including what leaves through portals.
struct RoomAcoustics {
    bool active = false;
    float volume = 0.f;
    float surfaceArea = 0.f;
    float meanFreePath = 0.f;
    Bands absorptionArea{};
    Bands decayRate{};
    Bands t60{};
};

// Energy exchange of the diffuse field: room `from` loses rate[band] * E to room `to` each second.
struct RoomLink {
    std::uint32_t from = kNone;
    std::uint32_t to = kNone;
    Bands rate{};
};

struct AcousticSnapshot {
    std::uint64_t tick = 0;
    Vec3 listener = vec3(0, 0, 0);
    std::uint32_t listenerRoom = kNone;
    std::uint32_t sourceCount = 0;
    std::uint32_t roomCount = 0;
    std::uint32_t linkCount = 0;
    SourceAcoustics sources[kMaxSources];
    RoomAcoustics rooms[kMaxRooms];
    RoomLink links[kMaxRoomLinks];
};

// One writer thread, one reader thread. Wait-free for both: each side only swaps an index with one atomic exchange.
template <class T>
class SnapshotExchange {
public:
    SnapshotExchange() : slots_(new Slots) {}

    T& writeSlot() { return slots_->values[back_]; }

    void publish() {
        const std::uint32_t previous = shared_.exchange(back_ | kFresh, std::memory_order_acq_rel);
        back_ = previous & kIndexMask;
    }

    // The newest published value, or the same one as last time if nothing new arrived. Valid until the next call to read.
    const T& read() {
        if (shared_.load(std::memory_order_acquire) & kFresh) front_ = shared_.exchange(front_, std::memory_order_acq_rel) & kIndexMask;
        return slots_->values[front_];
    }

private:
    static constexpr std::uint32_t kFresh = 4, kIndexMask = 3;
    struct Slots {
        T values[3];
    };
    std::unique_ptr<Slots> slots_;
    alignas(64) std::atomic<std::uint32_t> shared_{1};
    alignas(64) std::uint32_t back_ = 0;
    alignas(64) std::uint32_t front_ = 2;
};

}  // namespace mspatial::acoustics
