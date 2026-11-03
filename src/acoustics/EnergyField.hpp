// The listener's energy field: how much reflected energy reaches the listener, from which direction, at which moment, per octave band.
//
// Paths leave the listener, bounce through the live scene (specular or diffuse by each surface's scattering) and are scored by a receiver
// sphere around the listener: every segment that crosses it deposits energy * chord length into the bin of its arrival time and the
// direction it arrives from. By reciprocity this is the response of the room to a source at the listener, which is what the late field of
// any source in the same room converges to. The estimate is normalised to the free-field direct sound at 1 m, the unit the rest of the
// module uses, so it plugs into the same send levels as the feedback reverb.
//
// Tracing is streamed: each path advances one segment per batch, so the field costs a fixed number of rays per tick whatever its length,
// and the backend never waits. Old evidence fades with a time constant, so a moving door or a broken wall shows up within seconds.
#pragma once

#include "acoustics/AcousticScene.hpp"
#include "acoustics/AcousticTrace.hpp"

#include <array>
#include <span>
#include <vector>

namespace mspatial::acoustics {

// Arrival directions are binned onto the vertices of an icosahedron: twelve directions, each about 63 degrees from its neighbours.
inline constexpr std::uint32_t kFieldDirections = 12;
const std::array<Vec3, kFieldDirections>& fieldDirections();
std::uint32_t nearestFieldDirection(Vec3 direction);

struct EnergyField {
    std::uint32_t bins = 0;
    float binSeconds = 0.01f;
    std::uint32_t room = kNone;   // The listener's room when the evidence was gathered.
    double paths = 0.0;           // Effective number of paths behind the estimate, after fading.
    // Energy per bin relative to the direct sound at 1 m, laid out [direction][band][bin]. A bin holds energy, not power: summing the bins
    // gives the total reflected energy of an impulse.
    std::vector<float> energy;

    float& at(std::uint32_t direction, std::uint32_t band, std::uint32_t bin) { return energy[(direction * kBandCount + band) * bins + bin]; }
    float at(std::uint32_t direction, std::uint32_t band, std::uint32_t bin) const { return energy[(direction * kBandCount + band) * bins + bin]; }
    // All directions together.
    float band(std::uint32_t band, std::uint32_t bin) const;
    // Seconds for 60 dB of decay in each band from the backward-integrated field, fitted between -5 dB and as far as the field reaches
    // (at most -35 dB). Zero where the field holds too little to fit.
    Bands t60() const;
};

struct EnergyTracerConfig {
    std::uint32_t paths = 384;      // Paths in flight; each advances one segment per tick.
    float seconds = 1.5f;           // Length of the field.
    float binSeconds = 0.01f;
    float receiverRadius = 0.75f;   // Metres, at least. setRoomVolume widens it in large rooms.
    float timeConstant = 3.f;       // Seconds for old evidence to fade to 1/e.
    // Every path lives to the end of the field however weak it gets: the field's shape matters at every moment, and a path's cost is one ray
    // per tick whether it is loud or not. Only a path this far down in every band (-150 dB) ends early.
    float energyFloor = 1e-15f;
    float speedOfSound = kDefaultSpeedOfSound;
    Bands air = airAttenuation();
};

// The ray budget of one batch: paths whose segments ride in it.
struct EnergyJob {
    std::uint32_t path;
    std::uint32_t ray;
    std::uint32_t generation;
};

class EnergyTracer {
public:
    void configure(const EnergyTracerConfig& config);
    const EnergyTracerConfig& config() const { return config_; }
    bool enabled() const { return config_.paths > 0; }

    // The receiver grows with the room so a large hall gathers as much evidence per ray as a small one.
    void setRoomVolume(float volume);
    // A new room starts from nothing: evidence from another room describes another field.
    void setRoom(std::uint32_t room);

    // Queues at most budget rays, one per idle path, into the batch's closest-hit rays.
    std::uint32_t issue(TraceBatch& batch, std::uint32_t budget, Vec3 listener, float selfHitOffset, std::vector<EnergyJob>& jobs);
    void consume(const Scene& scene, const TraceBatch& batch, std::span<const EnergyJob> jobs, Vec3 listener);
    // Fades old evidence by dt seconds.
    void age(float dt);
    void reset();
    // With spawning off, a path that ends is not replaced, so the paths in flight run out: an offline bake drains them so that every path
    // it counts has reached the end of the field.
    void setSpawning(bool spawning) { spawning_ = spawning; }
    std::uint32_t livePaths() const;
    double emitted() const { return emitted_; }

    // The current estimate.
    void read(EnergyField& out) const;
    std::uint64_t segmentsTraced() const { return segments_; }

private:
    struct Path {
        Vec3 origin = vec3(0, 0, 0);
        Vec3 direction = vec3(0, 0, 1);
        float travelled = 0.f;
        Bands energy{};
        std::uint32_t order = 0;
        std::uint32_t generation = 0;
        bool live = false;
        bool inFlight = false;
    };

    void spawn(Path& path, Vec3 listener);

    EnergyTracerConfig config_;
    std::uint32_t bins_ = 0;
    float radius_ = 0.75f;
    std::uint32_t room_ = kNone;
    std::vector<Path> paths_;
    std::vector<double> sums_;  // Unnormalised chord * energy, [direction][band][bin].
    double emitted_ = 0.0;
    std::uint32_t generation_ = 1;
    std::uint32_t sequence_ = 0;
    std::uint32_t rng_ = 0x9e3779b9u;
    std::uint64_t segments_ = 0;
    bool spawning_ = true;
};

}  // namespace mspatial::acoustics
