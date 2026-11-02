#include "acoustics/RoomAcoustics.hpp"

#include "acoustics/AcousticScene.hpp"

#include <algorithm>
#include <numbers>

namespace mspatial::acoustics {

Bands airAttenuation(float celsius, float relativeHumidityPercent, float pressureKilopascal) {
    // ISO 9613-1 equations 2 to 5. The standard's result is in dB per metre of sound pressure level; m below is per metre of energy.
    const double T = double(celsius) + 273.15, T0 = 293.15, T01 = 273.16;
    const double pa = double(pressureKilopascal) / 101.325;
    const double saturation = std::pow(10.0, -6.8346 * std::pow(T01 / T, 1.261) + 4.6151);
    const double h = double(relativeHumidityPercent) * saturation / pa;
    const double frO = pa * (24.0 + 4.04e4 * h * (0.02 + h) / (0.391 + h));
    const double frN = pa * std::pow(T / T0, -0.5) * (9.0 + 280.0 * h * std::exp(-4.170 * (std::pow(T / T0, -1.0 / 3.0) - 1.0)));
    Bands out;
    for (std::size_t b = 0; b < kBandCount; ++b) {
        const double f = kBandCenterHz[b], f2 = f * f;
        const double dbPerMetre = 8.686 * f2 *
            (1.84e-11 / pa * std::sqrt(T / T0) +
             std::pow(T / T0, -2.5) * (0.01275 * std::exp(-2239.1 / T) / (frO + f2 / frO) + 0.1068 * std::exp(-3352.0 / T) / (frN + f2 / frN)));
        out[b] = float(dbPerMetre / (10.0 * std::log10(std::numbers::e)));
    }
    return out;
}

float equivalentAbsorption(DecayFormula formula, float surfaceArea, float absorptionArea) {
    if (formula == DecayFormula::Sabine || surfaceArea <= 0.f) return absorptionArea;
    const float alpha = std::clamp(absorptionArea / surfaceArea, 0.f, 0.999f);
    return -surfaceArea * std::log1p(-alpha);
}

float decayRate(float volume, float equivalentArea, float air, float speedOfSound) {
    return speedOfSound * (equivalentArea + 4.f * air * volume) / (4.f * volume);
}

float t60FromRate(float rate) { return rate > 0.f ? 6.f * std::numbers::ln10_v<float> / rate : INFINITY; }

ClosedMeasure measureClosedMesh(const Mesh& mesh, std::span<const Material> materials) {
    ClosedMeasure out;
    double volume = 0.0;
    for (std::size_t t = 0; t < mesh.triangleFacet.size(); ++t) {
        const Vec3 a = mesh.vertices[mesh.indices[t * 3]], b = mesh.vertices[mesh.indices[t * 3 + 1]], c = mesh.vertices[mesh.indices[t * 3 + 2]];
        volume += double(dot(a, cross(b, c))) / 6.0;
        const float area = 0.5f * length(cross(b - a, c - a));
        out.surfaceArea += area;
        const Material& material = materials[mesh.facets[mesh.triangleFacet[t]].material];
        for (std::size_t band = 0; band < kBandCount; ++band) out.absorptionArea[band] += area * material.absorption[band];
    }
    out.volume = float(std::fabs(volume));
    return out;
}

std::vector<float> CoupledRooms::steadyState(std::span<const float> powers) const {
    const std::size_t n = size();
    std::vector<double> a(n * n), x(powers.begin(), powers.end());
    for (std::size_t i = 0; i < n; ++i) {
        double diagonal = absorption[i];
        for (std::size_t j = 0; j < n; ++j) {
            diagonal += coupling[i * n + j];
            if (i != j) a[i * n + j] = -double(speedOfSound) / 4.0 * coupling[i * n + j];
        }
        a[i * n + i] = double(speedOfSound) / 4.0 * diagonal;
    }
    for (std::size_t k = 0; k < n; ++k) {  // Gaussian elimination; the matrix is diagonally dominant, so no pivoting is needed.
        for (std::size_t i = k + 1; i < n; ++i) {
            const double f = a[i * n + k] / a[k * n + k];
            for (std::size_t j = k; j < n; ++j) a[i * n + j] -= f * a[k * n + j];
            x[i] -= f * x[k];
        }
    }
    for (std::size_t k = n; k-- > 0;) {
        for (std::size_t j = k + 1; j < n; ++j) x[k] -= a[k * n + j] * x[j];
        x[k] /= a[k * n + k];
    }
    return {x.begin(), x.end()};
}

CoupledModes decayModes(const CoupledRooms& rooms) {
    const std::size_t n = rooms.size();
    std::vector<double> s(n * n, 0.0), q(n * n, 0.0);
    const double c4 = double(rooms.speedOfSound) / 4.0;
    for (std::size_t i = 0; i < n; ++i) {
        double diagonal = rooms.absorption[i];
        for (std::size_t j = 0; j < n; ++j) {
            diagonal += rooms.coupling[i * n + j];
            if (i != j) s[i * n + j] = -c4 * rooms.coupling[i * n + j] / std::sqrt(double(rooms.volumes[i]) * rooms.volumes[j]);
        }
        s[i * n + i] = c4 * diagonal / rooms.volumes[i];
        q[i * n + i] = 1.0;
    }
    for (int sweep = 0; sweep < 32; ++sweep) {  // Cyclic Jacobi.
        double off = 0.0;
        for (std::size_t p = 0; p < n; ++p)
            for (std::size_t r = p + 1; r < n; ++r) off += s[p * n + r] * s[p * n + r];
        if (off < 1e-20) break;
        for (std::size_t p = 0; p < n; ++p) {
            for (std::size_t r = p + 1; r < n; ++r) {
                const double apr = s[p * n + r];
                if (std::fabs(apr) < 1e-14) continue;
                const double theta = (s[r * n + r] - s[p * n + p]) / (2.0 * apr);
                const double t = std::copysign(1.0, theta) / (std::fabs(theta) + std::sqrt(theta * theta + 1.0));
                const double c = 1.0 / std::sqrt(t * t + 1.0), si = t * c;
                for (std::size_t k = 0; k < n; ++k) {
                    const double skp = s[k * n + p], skr = s[k * n + r];
                    s[k * n + p] = c * skp - si * skr;
                    s[k * n + r] = si * skp + c * skr;
                }
                for (std::size_t k = 0; k < n; ++k) {
                    const double spk = s[p * n + k], srk = s[r * n + k];
                    s[p * n + k] = c * spk - si * srk;
                    s[r * n + k] = si * spk + c * srk;
                }
                for (std::size_t k = 0; k < n; ++k) {
                    const double qkp = q[k * n + p], qkr = q[k * n + r];
                    q[k * n + p] = c * qkp - si * qkr;
                    q[k * n + r] = si * qkp + c * qkr;
                }
            }
        }
    }
    CoupledModes out;
    out.rates.resize(n);
    for (std::size_t i = 0; i < n; ++i) out.rates[i] = float(s[i * n + i]);
    out.shape = std::move(q);
    out.volumes = rooms.volumes;
    return out;
}

float CoupledModes::impulseEnergy(std::size_t from, std::size_t to, float t) const {
    // With x = sqrt(V) E the system is symmetric: x(t) = Q exp(-L t) Q^T x(0), and x(0) = Q_energy / sqrt(V_from).
    const std::size_t n = rates.size();
    double x = 0.0;
    for (std::size_t m = 0; m < n; ++m) x += shape[to * n + m] * shape[from * n + m] * std::exp(-double(rates[m]) * t);
    return float(x / std::sqrt(double(volumes[from]) * volumes[to]));
}

Bands apertureChainEnergy(std::span<const ApertureLeg> apertures, float lastDistance, float speedOfSound) {
    Bands out = uniformBands(0.f);
    if (apertures.empty()) return out;
    float unfolded = lastDistance;
    for (const ApertureLeg& a : apertures) unfolded += a.distance;
    const float freeField = 1.f / std::max(unfolded * unfolded, 1e-4f);
    const float pi = std::numbers::pi_v<float>;
    for (std::size_t band = 0; band < kBandCount; ++band) {
        const float k = 2.f * pi * kBandCenterHz[band] / speedOfSound;
        const auto lobe = [&](float cosine, float area) {
            const float ka = k * std::sqrt(area / pi);
            const float weight = ka / (1.f + ka);
            return (1.f - weight) + weight * cosine;
        };
        // Referenced to W / (4 pi): the first aperture collects A cos / d^2, and each later step spreads over pi d^2.
        const ApertureLeg& first = apertures[0];
        double energy = first.area * first.open[band] * lobe(first.cosIn, first.area) / std::max(double(first.distance) * first.distance, 1e-4);
        for (std::size_t i = 1; i < apertures.size(); ++i) {
            const ApertureLeg& previous = apertures[i - 1];
            const ApertureLeg& next = apertures[i];
            energy *= lobe(previous.cosOut, previous.area) * next.area * next.open[band] * lobe(next.cosIn, next.area) /
                      (pi * std::max(double(next.distance) * next.distance, 1e-4));
        }
        const ApertureLeg& last = apertures.back();
        energy *= lobe(last.cosOut, last.area) / (pi * std::max(double(lastDistance) * lastDistance, 1e-4));
        float limit = freeField;
        for (const ApertureLeg& a : apertures) limit *= a.open[band];
        out[band] = std::min(float(energy), limit);
    }
    return out;
}

}  // namespace mspatial::acoustics
