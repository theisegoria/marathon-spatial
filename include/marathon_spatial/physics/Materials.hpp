#pragma once
// Plain-data acoustic material presets. Copy one and change any field: nothing
// here is special to the library, and a host passes its own values in the same
// structs. The elastic figures are typical handbook values for the bulk
// material, and the texture figures are starting points for tuning. None of
// them has been measured from a sample.
#include "marathon_spatial/physics/Modal.hpp"
#include <span>
#include <string_view>

namespace mspatial::physics {
struct SurfaceTexture {
    double cornerWavelengthM{.001};    // SurfaceProfile corner: finer texture falls off
    double slideFluctuation{.25};      // force RMS / friction force, sliding
    double rollFluctuation{.06};       // force RMS / normal force, rolling
    double stickStiffnessNPerM{};      // tangential contact stiffness for stick-slip; 0 disables it
};
struct MaterialPreset {
    std::string_view name;
    ModalElastic elastic{};
    double restitution{.5};
    double staticFriction{.5},dynamicFriction{.4};
    SurfaceTexture texture{};
    double panelThicknessM{.01};       // when a thin surface of this material is modelled as a panel
};
[[nodiscard]] std::span<const MaterialPreset> materialPresets() noexcept;
[[nodiscard]] const MaterialPreset* findMaterialPreset(std::string_view name) noexcept;
[[nodiscard]] ModalMaterial modalMaterial(const MaterialPreset& preset);
}
