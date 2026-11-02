#include "marathon_spatial/physics/Materials.hpp"
#include <array>

namespace mspatial::physics {
namespace {
//  name          E Pa     rho    nu    loss     e    mu_s  mu_d  wave m   slide roll  stick N/m  panel m
constexpr std::array Presets{
    MaterialPreset{"steel",      {200e9,7850,.29,6e-4}, .60,.60,.45,{.0008,.18,.03,4e6},.012},
    MaterialPreset{"sheet-metal",{200e9,7850,.29,1.5e-3},.55,.60,.45,{.0006,.22,.05,2e6},.002},
    MaterialPreset{"aluminium",  {69e9,2700,.33,4e-4},  .55,.60,.45,{.0008,.20,.04,2e6},.006},
    MaterialPreset{"cast-iron",  {120e9,7200,.26,2e-3}, .50,.50,.40,{.0012,.25,.05,3e6},.02},
    MaterialPreset{"brass",      {100e9,8500,.34,1e-3}, .55,.50,.40,{.0006,.16,.03,2e6},.005},
    MaterialPreset{"copper",     {120e9,8960,.34,2e-3}, .50,.55,.42,{.0006,.16,.03,2e6},.004},
    MaterialPreset{"glass",      {70e9,2500,.22,1e-3},  .70,.50,.35,{.0002,.12,.02,1e6},.008},
    MaterialPreset{"ceramic",    {70e9,2400,.20,1.5e-3},.65,.50,.40,{.0005,.20,.04,1e6},.012},
    MaterialPreset{"concrete",   {30e9,2400,.20,2e-2},  .40,.70,.60,{.0015,.35,.12,0},.15},
    MaterialPreset{"brick",      {15e9,1900,.18,2e-2},  .40,.70,.60,{.002,.40,.14,0},.1},
    MaterialPreset{"granite",    {50e9,2700,.25,3e-3},  .50,.60,.50,{.0008,.25,.06,0},.03},
    MaterialPreset{"hardwood",   {12e9,720,.35,1e-2},   .50,.50,.35,{.0015,.25,.08,5e5},.02},
    MaterialPreset{"softwood",   {9e9,500,.35,1.5e-2},  .45,.50,.35,{.002,.28,.09,3e5},.02},
    MaterialPreset{"plywood",    {8e9,600,.30,1.5e-2},  .45,.45,.32,{.0015,.25,.08,3e5},.012},
    MaterialPreset{"acrylic",    {3.2e9,1190,.37,3e-2}, .50,.50,.40,{.0005,.20,.05,3e5},.004},
    MaterialPreset{"hdpe",       {1e9,950,.42,6e-2},    .40,.30,.25,{.0008,.25,.06,1e5},.006},
    MaterialPreset{"rubber",     {5e6,1100,.49,2e-1},   .80,1.0,.80,{.001,.30,.04,5e4},.01},
    MaterialPreset{"cardboard",  {3e9,700,.30,6e-2},    .20,.40,.30,{.001,.30,.06,0},.005},
};
}
std::span<const MaterialPreset> materialPresets() noexcept {return Presets;}
const MaterialPreset* findMaterialPreset(std::string_view name) noexcept {
    for(const auto& p:Presets)if(p.name==name)return &p;
    return nullptr;
}
ModalMaterial modalMaterial(const MaterialPreset& p){
    return modalMaterialFromLossFactor(p.elastic.densityKgPerM3,p.elastic.youngModulusPa,p.elastic.poissonRatio,p.elastic.lossFactor);
}
}
