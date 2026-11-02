#pragma once
// Contact-driven excitation: impacts, sliding (scrape, stick-slip) and rolling.
//
// A surface is a height field fixed in space, read at the contact position, so
// the same path over the same spot gives the same force. Sliding reads it at the
// slip speed. Rolling reads it through the Hertz contact patch, which is why a
// ball sounds smoother than a scrape over the same surface, and adds a
// once-per-turn component for out-of-round bodies. Forces feed ModalSynth as
// external drives, and modes respond linearly, so two contacts on one body are
// two voices that sum correctly.
//
// ContactAudio is single owner. Registration allocates. Everything else,
// including render, does not allocate and is deterministic: the same calls at the
// same frames give the same bits, for any render block size.
#include "marathon_spatial/physics/Materials.hpp"
#include "marathon_spatial/physics/Modal.hpp"
#include "marathon_spatial/physics/Source.hpp"
#include <cstdint>
#include <memory>
#include <optional>
#include <span>
#include <vector>

namespace mspatial::physics {
// Zero-mean, unit-RMS periodic height field with exponential autocorrelation and
// a filter pyramid. The height spectrum is flat at wavelengths above
// cornerWavelengthM and falls at 12 dB per octave below it, so a slide at v m/s
// has its corner at v / cornerWavelengthM Hz. 2^sizeLog2 points, spaced
// cornerWavelengthM / 32, before the surface repeats.
class SurfaceProfile {
public:
    SurfaceProfile(double cornerWavelengthM,std::uint32_t seed,std::uint32_t sizeLog2=18);
    // Height at xMeters, averaged over widthMeters (a contact patch or a sampling step).
    [[nodiscard]] float read(double xMeters,double widthMeters) const noexcept;
    [[nodiscard]] double cornerWavelengthM() const noexcept {return wavelength_;}
    [[nodiscard]] double spacingM() const noexcept {return spacing_;}
private:
    double wavelength_{},spacing_{};
    std::vector<std::vector<float>> levels_;
};

// J = (1 + e) m* v and the Hertz contact time for one impact.
struct ImpactInput {
    double normalSpeedMPerS{};
    double massAKg{},massBKg{};   // massBKg may be +infinity
    double radiusM{};             // effective contact radius
    double restitution{.5};
    ModalElastic a{},b{};
};
struct ImpactSolution {double impulseNs{},contactSeconds{};};
[[nodiscard]] ImpactSolution solveImpact(const ImpactInput& input);
// Hertz contact patch radius (3 N R / (4 E*))^(1/3).
[[nodiscard]] double hertzPatchRadiusM(double normalForceN,double radiusM,double effectiveModulusPa);
[[nodiscard]] double effectiveModulusPa(const ModalElastic& a,const ModalElastic& b);

enum class ContactKind : std::uint8_t {Slide,Roll,Source};
struct ContactForceParameters {
    ContactKind kind{ContactKind::Slide};
    double speedMPerS{};
    double normalForceN{};
    // Slide
    double dynamicFriction{.4},staticFriction{.5};
    double fluctuation{.25};            // slide: RMS / friction force. roll: RMS / normal force
    double stickStiffnessNPerM{};       // nonzero enables stick-slip below stribeckSpeedMPerS
    double stribeckSpeedMPerS{.2};
    // Roll
    double radiusM{.02};
    double effectiveModulusPa{1e10};
    double outOfRound{};                // once-per-flat force as a fraction of the normal force
    std::uint32_t flats{1};
    // Source: the host's block source, in newtons after forceScaleN. Not owned.
    BlockSource* source{};
    double forceScaleN{1};
    double smoothingSeconds{.02};       // speed and force follow their targets with this time constant
};
// Generates a contact force in newtons, zero mean. Usable alone.
class ContactForce {
public:
    ContactForce(double sampleRateHz,std::uint32_t seed);
    void set(const ContactForceParameters& parameters,const SurfaceProfile* surface) noexcept;
    void generate(std::span<float> force) noexcept;
    [[nodiscard]] double positionM() const noexcept {return position_;}
    [[nodiscard]] double speedMPerS() const noexcept {return speed_;}
    [[nodiscard]] double normalForceN() const noexcept {return normal_;}
private:
    ContactForceParameters p_{};
    const SurfaceProfile* surface_{};
    double fs_{},alpha_{},position_{},speed_{},normal_{},angle_{},stretch_{},stickMean_{},width_{};
    std::uint64_t frame_{};
    bool primed_{};
};

using SurfaceHandle=std::uint32_t;
struct ContactAudioSettings {
    ModalSynthSettings synth{};
    std::uint32_t maximumContacts{16};        // <= synth.maximumVoices
    std::uint32_t maximumScheduledHits{4096};
};
// One update of a sustained contact. contact is a nonzero key for this touch,
// not for the body: two touches on one body are two contacts. Speed, force and
// the model parameters are read each call. The force follows smoothly.
struct ContactUpdate {
    std::uint64_t contact{};
    ModalModelHandle model{};
    std::optional<std::uint32_t> sample;
    Vec3 localPosition{};
    std::optional<SurfaceHandle> surface;     // required for Slide and Roll
    ContactForceParameters force{};
    double rampSeconds{.02};
    double perceptualWeight{1};
    std::uint32_t maximumModes{};
};
struct ContactAudioStats {
    std::uint64_t impactsPlayed{},impactsRefused{},hitsDropped{},contactsRefused{};
};
class ContactAudio {
public:
    explicit ContactAudio(const ContactAudioSettings& settings);
    ~ContactAudio();
    ContactAudio(ContactAudio&&) noexcept;
    ContactAudio& operator=(ContactAudio&&) noexcept;
    [[nodiscard]] ModalModelHandle registerModel(std::shared_ptr<const ModalModel> model);
    [[nodiscard]] SurfaceHandle registerSurface(std::shared_ptr<const SurfaceProfile> surface);
    // Whether a call would be accepted, for a host checking before it crosses to the render thread.
    [[nodiscard]] bool valid(const ModalExcitation& excitation) const noexcept;
    [[nodiscard]] bool valid(const ContactUpdate& update) const noexcept;
    // Queues a hit delaySeconds from the next rendered frame, sample accurate.
    // Returns false when the queue is full or the hit is invalid. Hits at the same
    // frame play in call order.
    bool impact(const ModalExcitation& excitation,double delaySeconds=0) noexcept;
    // Starts or changes a sustained contact. Returns false if no contact slot is
    // free, and throws std::invalid_argument for an invalid update.
    bool update(const ContactUpdate& update);
    void end(std::uint64_t contact,double rampSeconds=.03) noexcept;
    void stop(std::uint64_t instance) noexcept;
    // Same layout as ModalSynth::render, for any frame count.
    void render(std::span<float> bus,std::span<float> voices,std::size_t frames);
    [[nodiscard]] std::uint32_t activeContacts() const noexcept;
    [[nodiscard]] std::uint64_t pendingHits() const noexcept;
    [[nodiscard]] const ContactForce* force(std::uint64_t contact) const noexcept;
    [[nodiscard]] const ContactAudioStats& stats() const noexcept;
    [[nodiscard]] ModalSynth& synth() noexcept;
    [[nodiscard]] const ModalSynth& synth() const noexcept;
    [[nodiscard]] double sampleRateHz() const noexcept;
private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};
}
