// SPDX-License-Identifier: Apache-2.0
#include "studio_bridge.h"

#include "../../../../include/marathon_spatial/physics/Modal.hpp"

#include <algorithm>
#include <exception>
#include <vector>

namespace physics = mspatial::physics;

extern "C" uint32_t studio_box_modes(const ms_physics_material* m, ms_vec3 size, float wall, double* hz, double* damping, uint32_t capacity) {
    if (!m || !hz || !damping) return 0;
    try {
        const physics::ModalElastic elastic{m->young_modulus_pa, m->density_kg_m3, m->poisson_ratio, m->loss_factor};
        const auto model = physics::analyticModalModel(physics::modalBodyForBox(physics::Vec3(size.x, size.y, size.z), wall, elastic));
        const auto count = static_cast<uint32_t>(std::min<std::size_t>(model.modes.size(), capacity));
        for (uint32_t i = 0; i < count; ++i) {
            hz[i] = model.modes[i].frequencyHz;
            damping[i] = model.modes[i].dampingRatio;
        }
        return count;
    } catch (const std::exception&) {
        return 0;
    }
}

extern "C" ms_acoustics* studio_acoustics_create(const studio_room_desc* room, const char* metallib, float sample_rate) {
    if (!room || !room->positions || !room->indices || !room->triangle_materials || !room->absorption || !room->scattering || !room->transmission ||
        !room->absorption_area || room->material_count == 0)
        return nullptr;
    try {
        ms_acoustics_config config;
        ms_acoustics_config_default(&config);
        config.metallib_path = metallib;
        config.sample_rate = sample_rate;
        ms_acoustics* acoustics = ms_acoustics_create(&config);
        if (!acoustics) return nullptr;

        std::vector<uint32_t> ids(room->material_count);
        for (uint32_t m = 0; m < room->material_count; ++m) {
            ms_acoustic_material material{};
            for (uint32_t band = 0; band < MS_ACOUSTIC_BANDS; ++band) {
                material.absorption[band] = room->absorption[m * MS_ACOUSTIC_BANDS + band];
                material.scattering[band] = room->scattering[m];
                material.transmission[band] = room->transmission[m];
            }
            ids[m] = ms_acoustics_material_add(acoustics, &material);
        }
        std::vector<uint32_t> triangles(room->index_count / 3);
        for (std::size_t t = 0; t < triangles.size(); ++t) triangles[t] = ids[std::min(room->triangle_materials[t], room->material_count - 1)];
        const uint32_t mesh = ms_acoustics_mesh_add(acoustics, room->positions, room->vertex_count, room->indices, room->index_count, triangles.data());
        const ms_acoustic_transform identity{{{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}}};
        (void)ms_acoustics_instance_add(acoustics, mesh, &identity);

        const ms_acoustic_box bounds{room->bounds_min, room->bounds_max};
        ms_acoustic_room desc{};
        desc.bounds = &bounds;
        desc.bound_count = 1;
        desc.probe = room->probe;
        desc.volume = room->volume;
        desc.surface_area = room->surface_area;
        std::copy_n(room->absorption_area, MS_ACOUSTIC_BANDS, desc.absorption_area);
        (void)ms_acoustics_room_add(acoustics, &desc);

        // Reverb delays come from a published room, so the first tick precedes the prepare.
        ms_acoustics_tick(acoustics, 1.f / 60.f);
        if (ms_acoustics_reverb_prepare(acoustics) != MS_OK) {
            ms_acoustics_destroy(acoustics);
            return nullptr;
        }
        return acoustics;
    } catch (const std::exception&) {
        return nullptr;
    }
}

struct studio_reverb {
    ms_acoustics* acoustics;
    uint32_t capacity;
    std::vector<float> store;
};

extern "C" studio_reverb* studio_reverb_create(ms_acoustics* acoustics, uint32_t max_frames) {
    if (!acoustics || max_frames == 0) return nullptr;
    try {
        return new studio_reverb{acoustics, max_frames, std::vector<float>(std::size_t(max_frames) * 3, 0.f)};
    } catch (const std::exception&) {
        return nullptr;
    }
}

extern "C" void studio_reverb_destroy(studio_reverb* reverb) { delete reverb; }

extern "C" void studio_reverb_bus(void* user, float* left, float* right, uint32_t frames, uint64_t) {
    auto& reverb = *static_cast<studio_reverb*>(user);
    if (frames > reverb.capacity) {
        std::fill_n(left, frames, 0.f);
        std::fill_n(right, frames, 0.f);
        return;
    }
    // A send carries the same mono on both sides.
    float* in = reverb.store.data();
    float* outLeft = in + reverb.capacity;
    float* outRight = outLeft + reverb.capacity;
    std::copy_n(left, frames, in);
    const float* inputs[1] = {in};
    float* lefts[1] = {outLeft};
    float* rights[1] = {outRight};
    ms_acoustics_reverb_process(reverb.acoustics, inputs, lefts, rights, frames);
    std::copy_n(outLeft, frames, left);
    std::copy_n(outRight, frames, right);
}
