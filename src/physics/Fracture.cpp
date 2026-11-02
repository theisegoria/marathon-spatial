#include "marathon_spatial/physics/Fracture.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <string>

namespace mspatial::physics {
namespace {
struct Rng {
    std::uint64_t state;
    double uniform() noexcept {state+=0x9e3779b97f4a7c15ull;std::uint64_t z=state;z=(z^(z>>30))*0xbf58476d1ce4e5b9ull;z=(z^(z>>27))*0x94d049bb133111ebull;z^=z>>31;return (double(z>>11)+.5)/9007199254740992.;}
    std::uint32_t below(std::uint32_t n) noexcept {return static_cast<std::uint32_t>(uniform()*n);}
};
[[noreturn]] void bad(const std::string& what){throw std::invalid_argument("fracture: "+what);}
}
DebrisBank makeDebrisBank(ContactAudio& audio,const MaterialPreset& material,double lo,double hi,std::uint32_t count,std::uint32_t modes){
    if(!std::isfinite(lo)||!(lo>0)||!std::isfinite(hi)||hi<lo)bad("bank mass range invalid");
    if(count<1||count>32)bad("bank size outside 1..32");
    DebrisBank bank;bank.material=material.elastic;bank.restitution=material.restitution;
    for(std::uint32_t i=0;i<count;++i){
        const double mass=count>1?lo*std::pow(hi/lo,double(i)/(count-1)):lo;
        // A chunk too small to ring in the audible band still ticks. Hold the model at the smallest
        // size whose lowest modes stay under 14 kHz. The hit keeps the fragment's real mass.
        const double shear=std::sqrt(material.elastic.youngModulusPa/(2*(1+material.elastic.poissonRatio)*material.elastic.densityKgPerM3));
        const double edge=std::max(std::cbrt(mass/material.elastic.densityKgPerM3),shear/(2*14000));
        // Rough chunks: not quite cubes, so the entries do not share one mode pattern.
        ModalBody body;body.shape=ModalShape::Block;body.size={float(edge*1.3),float(edge),float(edge/1.3)};body.material=material.elastic;
        body.maximumFrequencyHz=22000;body.minimumFrequencyHz=40;body.maximumModes=modes;body.samplesPerAxis=3;
        auto model=analyticModalModel(body);
        const auto samples=static_cast<std::uint32_t>(model.samples.size());
        bank.entries.push_back({audio.registerModel(std::make_shared<const ModalModel>(std::move(model))),mass,edge/2,samples});
    }
    return bank;
}
FracturePlan planFracture(const FractureEvent& e,const DebrisBank& bank){
    if(e.id==0)bad("id must be nonzero");
    if(!std::isfinite(e.parentMassKg)||!(e.parentMassKg>0))bad("parent mass must be positive");
    if(!std::isfinite(e.energyJ)||e.energyJ<0)bad("energy must be non-negative");
    if(!std::isfinite(e.kineticFraction)||e.kineticFraction<0||e.kineticFraction>1)bad("kinetic fraction outside [0, 1]");
    if(e.fragments<1||e.fragments>4096)bad("fragment count outside 1..4096");
    if(!std::isfinite(e.massExponent)||!(e.massExponent>0))bad("mass exponent must be positive");
    if(!std::isfinite(e.heightM)||e.heightM<0||!std::isfinite(e.gravityMPerS2)||!(e.gravityMPerS2>0))bad("height or gravity invalid");
    if(bank.entries.empty())bad("empty debris bank");
    if(!(e.floor.youngModulusPa>0))bad("floor material required");
    Rng rng{(std::uint64_t(e.seed)<<32)^(e.id*0x9e3779b97f4a7c15ull)};
    FracturePlan plan;
    // Stratified inverse-CDF of N(>m) ~ m^-b truncated to a 1:10^4 mass range, scaled to the parent's mass.
    const double b=e.massExponent,ratio=1e4;
    double total=0;
    for(std::uint32_t i=0;i<e.fragments;++i){
        const double u=(i+rng.uniform())/e.fragments;
        const double m=std::pow(1-u*(1-std::pow(ratio,-b)),-1/b);
        plan.fragmentMassKg.push_back(m);total+=m;
    }
    for(auto& m:plan.fragmentMassKg)m*=e.parentMassKg/total;
    // Speed: one shared specific energy with scatter, then scaled so the kinetic energy is exact.
    double energy=0;
    for(std::uint32_t i=0;i<e.fragments;++i){
        const double v=std::exp(.5*(rng.uniform()+rng.uniform()-1));
        plan.fragmentSpeedMPerS.push_back(v);energy+=.5*plan.fragmentMassKg[i]*v*v;
    }
    const double target=e.energyJ*e.kineticFraction;
    const double scale=energy>0?std::sqrt(target/energy):0;
    for(auto& v:plan.fragmentSpeedMPerS)v*=scale;
    plan.kineticEnergyJ=target;
    for(std::uint32_t i=0;i<e.fragments;++i){
        const double mass=plan.fragmentMassKg[i];
        std::uint32_t entry=0;double nearest=1e300;
        for(std::uint32_t k=0;k<bank.entries.size();++k){
            const double d=std::fabs(std::log(mass/bank.entries[k].massKg));
            if(d<nearest){nearest=d;entry=k;}
        }
        const double vertical=plan.fragmentSpeedMPerS[i]*(2*rng.uniform()-1); // isotropic: uniform vertical cosine
        const double g=e.gravityMPerS2;
        double speed=std::sqrt(vertical*vertical+2*g*e.heightM);
        double t=(vertical+speed)/g;
        const double radius=bank.entries[entry].radiusM*std::cbrt(mass/bank.entries[entry].massKg);
        for(std::uint32_t bounce=0;bounce<=e.maximumBounces&&speed>=e.minimumSpeedMPerS;++bounce){
            const auto solved=solveImpact({speed,mass,std::numeric_limits<double>::infinity(),radius,bank.restitution,bank.material,e.floor});
            plan.hits.push_back({t,i,entry,rng.below(std::max(bank.entries[entry].samples,1u)),speed,solved.impulseNs,solved.contactSeconds});
            speed*=bank.restitution;
            t+=2*speed/g;
        }
    }
    std::stable_sort(plan.hits.begin(),plan.hits.end(),[](const DebrisHit& a,const DebrisHit& c){return a.timeSeconds<c.timeSeconds;});
    return plan;
}
ModalExcitation debrisExcitation(const FractureEvent& e,const DebrisHit& h,const DebrisBank& bank,double weight) noexcept {
    ModalExcitation x;
    x.instance=((e.id*0x9e3779b97f4a7c15ull)&~std::uint64_t(63))|(h.bankEntry+1); // a bank has at most 32 entries
    x.model=bank.entries[h.bankEntry].model;x.sample=h.sample;x.impulseNs=h.impulseNs;
    x.contactSeconds=std::clamp(h.contactSeconds,MinimumContactSeconds,MaximumContactSeconds);x.perceptualWeight=weight;
    return x;
}
std::size_t scheduleFracture(ContactAudio& audio,const FractureEvent& e,const FracturePlan& plan,const DebrisBank& bank,double delay,double weight) noexcept {
    std::size_t accepted=0;
    for(const auto& h:plan.hits)accepted+=audio.impact(debrisExcitation(e,h,bank,weight),delay+h.timeSeconds)?1:0;
    return accepted;
}
}
