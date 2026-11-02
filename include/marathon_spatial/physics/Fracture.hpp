#pragma once
// Fracture and debris clusters. A break becomes a planned rain of hits: fragment
// masses from a truncated power law, speeds from the released kinetic energy,
// ballistic flight, then bounces that lose speed by the restitution. The hits go
// through a small bank of fragment models, one voice per bank entry. Resonators
// are linear, so hundreds of hits on one voice cost no more than a few.
#include "marathon_spatial/physics/Contact.hpp"
#include <cstdint>
#include <span>
#include <vector>

namespace mspatial::physics {
struct DebrisBank {
    struct Entry {ModalModelHandle model{};double massKg{},radiusM{};std::uint32_t samples{};};
    std::vector<Entry> entries;   // ascending mass
    ModalElastic material{};
    double restitution{.5};
};
// count models, log-spaced from minimumMassKg to maximumMassKg, registered with audio.
[[nodiscard]] DebrisBank makeDebrisBank(ContactAudio& audio,const MaterialPreset& material,double minimumMassKg,double maximumMassKg,std::uint32_t count,std::uint32_t modesPerEntry=12);

struct FractureEvent {
    std::uint64_t id{};               // nonzero
    double parentMassKg{};
    double energyJ{};                 // released by the break
    double kineticFraction{.4};       // share of energyJ that becomes fragment motion
    std::uint32_t fragments{};
    double massExponent{.8};          // b of N(>m) ~ m^-b
    double heightM{1};                // fragment start above the floor
    double gravityMPerS2{9.81};
    std::uint32_t maximumBounces{5};
    double minimumSpeedMPerS{.05};    // a bounce slower than this ends the fragment
    ModalElastic floor{};             // the surface it lands on
    std::uint32_t seed{1};
};
struct DebrisHit {
    double timeSeconds{};
    std::uint32_t fragment{},bankEntry{},sample{};
    double speedMPerS{},impulseNs{},contactSeconds{};
};
struct FracturePlan {
    std::vector<double> fragmentMassKg,fragmentSpeedMPerS;
    std::vector<DebrisHit> hits;      // ascending time
    double kineticEnergyJ{};
};
// Deterministic for one event and bank.
[[nodiscard]] FracturePlan planFracture(const FractureEvent& event,const DebrisBank& bank);
// One planned hit as an excitation: every hit on a bank entry shares one instance.
[[nodiscard]] ModalExcitation debrisExcitation(const FractureEvent& event,const DebrisHit& hit,const DebrisBank& bank,double perceptualWeight=1) noexcept;
// Queues the plan delaySeconds from now. Returns the hits accepted.
std::size_t scheduleFracture(ContactAudio& audio,const FractureEvent& event,const FracturePlan& plan,const DebrisBank& bank,double delaySeconds=0,double perceptualWeight=1) noexcept;
}
