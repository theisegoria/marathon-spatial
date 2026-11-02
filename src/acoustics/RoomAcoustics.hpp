// Closed-form room acoustics. Everything here is analytic and independent of the ray backend.
#pragma once

#include "acoustics/AcousticTypes.hpp"

namespace mspatial::acoustics {

enum class DecayFormula { Sabine, Eyring };

// Equivalent absorption area in m^2: the area of open window with the same effect. Sabine sums alpha S; Eyring uses -S ln(1 - mean alpha),
// which stays right for dead rooms where Sabine overestimates the decay time.
float equivalentAbsorption(DecayFormula formula, float surfaceArea, float absorptionArea);

// Energy decay rate per second of the diffuse field: c (A + 4 m V) / (4 V).
float decayRate(float volume, float equivalentArea, float airAttenuation, float speedOfSound);

// Time for energy to fall 60 dB at a given energy decay rate.
float t60FromRate(float rate);

// Volume, surface and absorption area of a closed triangle mesh, by the divergence theorem. Winding must be consistent.
struct ClosedMeasure {
    float volume = 0.f;
    float surfaceArea = 0.f;
    Bands absorptionArea{};
};
struct Mesh;
struct Material;
ClosedMeasure measureClosedMesh(const Mesh& mesh, std::span<const Material> materials);

// Two-room and many-room diffuse field. The energy densities E follow
//   V_i dE_i/dt = -(c/4) (A_i + sum_j S_ij) E_i + (c/4) sum_j S_ij E_j + P_i
// with A_i the absorption of room i and S_ij the coupling area. Volume-weighting makes the system symmetric, so its modes are real.
struct CoupledRooms {
    std::vector<float> volumes;
    std::vector<float> absorption;   // A_i in m^2.
    std::vector<float> coupling;     // S_ij in m^2, n by n, symmetric, zero diagonal.
    float speedOfSound = kDefaultSpeedOfSound;

    std::size_t size() const { return volumes.size(); }
    // Steady energy densities for constant powers.
    std::vector<float> steadyState(std::span<const float> powers) const;
};

struct CoupledModes {
    std::vector<float> rates;   // Energy decay rate of each mode, per second.
    std::vector<double> shape;  // n by n: symmetric eigenvectors in columns.
    std::vector<float> volumes;

    // Energy density in room `to` at time t after an impulse of energy Q = 1 in room `from`.
    float impulseEnergy(std::size_t from, std::size_t to, float t) const;
};
CoupledModes decayModes(const CoupledRooms& rooms);

// Energy gain of a point source heard through a chain of small apertures, relative to the free-field level at 1 m. The source floods
// the first aperture, each aperture re-radiates what it received as a Lambert radiator, and low bands lose the cosine law because
// an aperture small against the wavelength radiates evenly. The result never exceeds the free-field level at the unfolded distance.
struct ApertureLeg {
    float distance = 0.f;  // Length of the leg arriving at this aperture from the previous point.
    float cosIn = 1.f;     // Incidence at this aperture.
    float cosOut = 1.f;    // Exit angle toward the next point.
    float area = 0.f;
    Bands open = uniformBands(1.f);  // Energy fraction let through in each band.
};
Bands apertureChainEnergy(std::span<const ApertureLeg> apertures, float lastDistance, float speedOfSound);

}  // namespace mspatial::acoustics
