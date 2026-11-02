#pragma once
// Physics-driven modal sound.
//
// Cook time: a closed triangle mesh (or a convex hull given as its triangles)
// is filled with a regular voxel lattice split into five P1 tetrahedra per
// cell, linear-elastic stiffness and lumped mass are assembled densely, and
// the generalized eigenproblem K u = lambda M u is solved through Accelerate
// LAPACK (mass-scaled to a standard problem, dsyevr). The six rigid modes are
// verified and discarded; elastic modes inside the audible band become a
// versioned `ModalModel` whose canonical encoding is named by SHA-256.
// `analyticModalModel` gives closed-form bars, plates and blocks without a mesh.
//
// Run time: `ModalSynth` binds voices to sounding instances of registered
// models. An excitation (contact sample or local position, normal impulse in
// N s, half-sine contact duration) sets per-mode damped sinusoids rendered by
// two-pole recursive resonators with double state. A drive feeds the same
// resonators a continuous contact force, for scrapes and rolling.
//
// Nothing here opens an audio device. Buffers go to whoever calls render.
#include <array>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <span>
#include <string_view>
#include <vector>

namespace mspatial::physics {
// Plain float vector. Converts both ways with any host vector that has float
// x, y, z members, so a game passes its own type.
template<class V> concept XyzVector=requires(const V& v){{v.x}->std::convertible_to<float>;{v.y}->std::convertible_to<float>;{v.z}->std::convertible_to<float>;};
struct Vec3 {
    float x{},y{},z{};
    constexpr Vec3()=default;
    constexpr Vec3(float a,float b,float c):x(a),y(b),z(c){}
    template<XyzVector V> requires(!std::same_as<V,Vec3>) constexpr Vec3(const V& v):x(v.x),y(v.y),z(v.z){}
    template<XyzVector V> requires(!std::same_as<V,Vec3>&&std::constructible_from<V,float,float,float>) constexpr operator V() const {return V(x,y,z);}
    constexpr Vec3 operator+(Vec3 b) const {return {x+b.x,y+b.y,z+b.z};}
    constexpr Vec3 operator-(Vec3 b) const {return {x-b.x,y-b.y,z-b.z};}
};

// Linear-elastic material and Rayleigh damping C = alpha M + beta K, giving the
// modal damping ratio xi = (alpha / omega + beta omega) / 2.
struct ModalMaterial {
    double densityKgPerM3{};
    double youngModulusPa{};
    double poissonRatio{};           // [0, 0.495]; P1 tetrahedra lock (overestimate frequency) as nu approaches 0.5
    double rayleighAlphaPerSecond{}; // mass-proportional damping (>= 0)
    double rayleighBetaSeconds{};    // stiffness-proportional damping (>= 0)
};
// Throws std::invalid_argument, naming `who`, for values the analysis rejects.
void validateModalMaterial(const ModalMaterial& material,std::string_view who);
// Rayleigh terms that give xi = eta / 2 at 200 Hz and at 5 kHz.
[[nodiscard]] ModalMaterial modalMaterialFromLossFactor(double densityKgPerM3,double youngModulusPa,double poissonRatio,double lossFactor);

// Closed, consistently outward-wound triangle mesh in model-local metres.
struct ModalMesh {
    std::span<const Vec3> vertices;
    std::span<const std::uint32_t> indices; // 3 per triangle
};
struct ModalAnalysisSettings {
    static constexpr std::uint32_t MaximumDegreesOfFreedom=3600; // dense solve bound (~104 MB matrix)
    static constexpr std::uint32_t MaximumModes=256;
    static constexpr std::uint32_t MaximumSurfaceSamples=1024;
    double latticeSpacingMeters{};   // cubic cell edge; required
    double minimumFrequencyHz{20},maximumFrequencyHz{20000};
    std::uint32_t maximumModes{64};
    std::uint32_t maximumSurfaceSamples{128};
    std::uint32_t maximumDegreesOfFreedom{MaximumDegreesOfFreedom};
};

struct ModalMode {
    double angularFrequency{}; // undamped, rad/s
    double frequencyHz{};      // angularFrequency / 2 pi
    double dampingRatio{};     // [0, 1)
    double surfaceRms{};       // RMS normal mass-normalized shape over the surface samples (kg^-1/2)
};
struct ModalSurfaceSample {
    Vec3 position{}; // closest point on the source mesh to the lattice node
    Vec3 normal{};   // outward unit normal of the nearest source triangles
    std::uint32_t node{};
};
struct ModalModel {
    static constexpr std::uint32_t FormatVersion=1;
    ModalMaterial material{};
    double latticeSpacingMeters{},minimumFrequencyHz{},maximumFrequencyHz{};
    std::uint32_t latticeNodes{},latticeTetrahedra{};
    double latticeMassKg{};
    Vec3 boundsMinimum{},boundsMaximum{}; // source mesh bounds
    std::array<unsigned char,32> sourceMeshSHA256{};
    std::vector<ModalMode> modes;              // ascending frequency
    std::vector<ModalSurfaceSample> samples;
    // samples.size() x modes.size(), sample-major: mass-normalized mode shape at
    // the sample node projected on the sample normal (kg^-1/2).
    std::vector<float> gains;
    std::array<unsigned char,32> identity{};   // SHA-256 of the canonical encoding
    [[nodiscard]] float gain(std::size_t sample,std::size_t mode) const {return gains[sample*modes.size()+mode];}
};
struct ModalAnalysisReport {
    std::uint64_t gridCells{};
    std::uint32_t latticeCells{},latticeNodes{},latticeTetrahedra{},degreesOfFreedom{},surfaceNodes{};
    double meshVolume{},latticeVolume{},latticeMassKg{};
    std::array<double,6> rigidEigenvalues{}; // (rad/s)^2, expected ~0
    double rigidTolerance{};                  // admitted |lambda| for a rigid mode
    double rigidSubspaceResidual{};           // max distance of the six vectors from the analytic rigid subspace
    double firstElasticEigenvalue{};
    std::uint32_t eigenpairsComputed{},solverPasses{};
    std::uint32_t modesBelowBand{},overdampedModes{};
    double latticeSeconds{},solveSeconds{}; // wall clock (voxelization; assembly + eigensolve); diagnostics only
};
struct ModalAnalysisResult {ModalModel model;ModalAnalysisReport report;};
// Throws std::invalid_argument for an open/inverted/degenerate mesh, invalid
// material or settings, a lattice that is empty, disconnected or above the DOF
// bound (the message gives a coarser spacing), or a rigid-mode check failure.
[[nodiscard]] ModalAnalysisResult analyzeModal(const ModalMesh& mesh,const ModalMaterial& material,const ModalAnalysisSettings& settings);

// Canonical little-endian encoding with a trailing SHA-256 of all prior bytes.
// encode validates first; decode rejects any malformed, non-canonical or
// tampered byte string. identity is computed, never trusted from the input.
void validateModalModel(const ModalModel& model);
[[nodiscard]] std::vector<unsigned char> encodeModalModel(const ModalModel& model);
[[nodiscard]] ModalModel decodeModalModel(std::span<const unsigned char> bytes);
[[nodiscard]] std::array<unsigned char,32> modalModelIdentity(const ModalModel& model);
// Nearest surface sample (Euclidean, lowest index on ties) to a model-local point.
[[nodiscard]] std::uint32_t nearestModalSample(const ModalModel& model,Vec3 localPosition);

// Closed-form model of a box-like body, with no mesh and no eigensolve: a free
// bar (Euler-Bernoulli, in both bending planes), a simply supported rectangular
// Kirchhoff plate, or a block (shear box modes). Plates and bars are exact for
// their idealization. Blocks are an estimate, and analyzeModal is the reference
// when one matters.
enum class ModalShape : std::uint8_t {Bar,Plate,Block};
struct ModalElastic {
    double youngModulusPa{},densityKgPerM3{},poissonRatio{},lossFactor{};
};
struct ModalBody {
    ModalShape shape{ModalShape::Plate};
    Vec3 size{};                 // longest, middle, thickness (m), centred on the origin
    ModalElastic material{};
    double airDampingPerSecond{0}; // adds xi = air / (2 omega)
    double minimumFrequencyHz{20},maximumFrequencyHz{20000};
    std::uint32_t maximumModes{32};
    std::uint32_t samplesPerAxis{4}; // [2, 12]
};
[[nodiscard]] ModalModel analyticModalModel(const ModalBody& body);
// Bar, plate or block by aspect ratio. wallThicknessM > 0 makes a thin-walled panel.
[[nodiscard]] ModalBody modalBodyForBox(Vec3 size,float wallThicknessM,const ModalElastic& material);

// Normalized spectrum of a unit-impulse half-sine force of the given duration,
// cos(omega tau / 2) / (1 - x^2) with x = omega tau / pi: 1 at DC, pi/4 at
// x = 1, zeros at odd x >= 3, then a 1/x^2 roll-off. Signed.
[[nodiscard]] double halfSinePulseSpectrum(double angularFrequency,double contactSeconds);
inline constexpr double MinimumContactSeconds=1e-6,MaximumContactSeconds=0.1;
// Hertz sphere-on-body impact duration estimate
// tau = 2.87 (m*^2 / (R* E*^2 v))^(1/5), with 1/E* = (1-nu_a^2)/E_a + (1-nu_b^2)/E_b.
// Pass infinite partner mass for a static body. Clamped to the contact bounds.
struct HertzContact {
    double relativeNormalSpeed{};   // m/s
    double massAKg{},massBKg{};     // massBKg may be +infinity
    double effectiveRadiusMeters{}; // contact curvature radius
    double youngModulusAPa{},poissonRatioA{},youngModulusBPa{},poissonRatioB{};
};
[[nodiscard]] double hertzContactSeconds(const HertzContact& contact);

struct ModalSynthSettings {
    static constexpr std::uint32_t MaximumVoices=64,MaximumModesPerVoice=256,MaximumBlockFrames=8192;
    double sampleRateHz{48000};
    std::uint32_t maximumVoices{16};
    std::uint32_t maximumModesPerVoice{64};
    std::uint32_t maximumBlockFrames{1024};
    // Output amplitude (outputGain-weighted surface RMS normal velocity, m/s)
    // under which a mode is culled; 0 disables culling.
    double cullAmplitude{1e-6};
    double outputGain{1};
};
using ModalModelHandle=std::uint32_t;
struct ModalExcitation {
    std::uint64_t instance{};               // nonzero persistent sounding-object key
    ModalModelHandle model{};
    std::optional<std::uint32_t> sample;    // contact sample; otherwise nearest to localPosition
    Vec3 localPosition{};
    double impulseNs{};                     // normal impulse magnitude, >= 0
    double contactSeconds{1e-3};            // half-sine force duration
    double perceptualWeight{1};             // > 0: scales energy when voices compete (distance gain at the listener)
    std::uint32_t maximumModes{};           // nonzero: only this many strongest modes (distance LOD)
};
// A continuous contact force on one surface sample. By default a zero-mean noise
// of forceRmsN whose spectrum falls above cutoffHz (surface texture read at the
// slip or rolling speed). With externalForce the host supplies the force per
// frame through setDriveForce, and forceRmsN only ranks voices. Repeat to change
// it; forceRmsN 0 releases the voice, which then rings down. The force ramps over
// rampSeconds between calls.
struct ModalDrive {
    std::uint64_t instance{};               // nonzero persistent sounding-object key
    ModalModelHandle model{};
    std::optional<std::uint32_t> sample;
    Vec3 localPosition{};
    double forceRmsN{};                     // >= 0
    double cutoffHz{1000};                  // > 0
    double rampSeconds{.02};                // >= 0
    double perceptualWeight{1};             // > 0
    std::uint32_t maximumModes{};           // nonzero: only this many most responsive modes
    bool externalForce{};
};
enum class ModalExciteOutcome : std::uint8_t {NewVoice,Merged,Stole,RejectedQuieter,Inaudible};
struct ModalExciteResult {
    ModalExciteOutcome outcome{ModalExciteOutcome::Inaudible};
    std::uint32_t voice{};          // slot, valid unless RejectedQuieter/Inaudible
    std::uint32_t sample{};         // resolved contact sample
    std::uint32_t modesExcited{};
    double excitationEnergy{};      // sum of squared output amplitudes / 2 of the accepted modes
    std::uint64_t stolenInstance{}; // when Stole
};
struct ModalVoiceState {
    bool active{};
    std::uint64_t instance{};
    ModalModelHandle model{};
    std::uint32_t activeModes{};
    double energy{};
    bool driven{};
};
// Single owner: not thread-safe. registerModel allocates (load time); excite,
// drive, stop and render do not allocate. Deterministic: identical call
// sequences give identical output bits.
class ModalSynth {
public:
    explicit ModalSynth(const ModalSynthSettings& settings);
    ~ModalSynth();
    ModalSynth(ModalSynth&&) noexcept;
    ModalSynth& operator=(ModalSynth&&) noexcept;
    [[nodiscard]] ModalModelHandle registerModel(std::shared_ptr<const ModalModel> model);
    // Starts at the next rendered frame. Merges into the instance's voice,
    // otherwise takes a free voice, otherwise steals the lowest-energy voice if
    // the new excitation carries more energy (each scaled by its perceptual
    // weight). A driven voice is not stolen by an impact.
    ModalExciteResult excite(const ModalExcitation& excitation);
    // Starts, changes or (forceRmsN 0) releases the instance's continuous force.
    // Takes a voice as excite does; only an idle or decaying voice can be stolen.
    ModalExciteResult drive(const ModalDrive& drive);
    // The force for the next render of an externalForce voice: at least that
    // many frames, read once and not kept. With none set the force is 0.
    void setDriveForce(std::uint32_t slot,const float* force) noexcept;
    [[nodiscard]] bool driven(std::uint32_t slot) const noexcept;
    [[nodiscard]] std::uint64_t instanceIn(std::uint32_t slot) const noexcept; // 0 when the slot is idle
    [[nodiscard]] std::uint32_t modelCount() const noexcept;
    [[nodiscard]] std::uint32_t sampleCount(ModalModelHandle model) const noexcept; // 0 for an unknown handle
    void stop(std::uint64_t instance) noexcept;
    // bus: frames samples (sum of voices). voices: empty, or maximumVoices x
    // frames, slot-major mono output per voice (inactive slots are zeroed).
    void render(std::span<float> bus,std::span<float> voices,std::size_t frames);
    [[nodiscard]] ModalVoiceState voice(std::uint32_t slot) const;
    [[nodiscard]] std::uint32_t activeVoices() const noexcept;
    [[nodiscard]] const ModalSynthSettings& settings() const noexcept;
private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};
}
