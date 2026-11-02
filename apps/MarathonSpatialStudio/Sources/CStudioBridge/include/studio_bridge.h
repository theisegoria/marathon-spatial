/* SPDX-License-Identifier: Apache-2.0
 *
 * The library's C API for Swift, plus what Swift should not do itself: the modes of a box for display, one room's acoustic
 * scene built from flat arrays, and the render-thread side of its late reverb.
 */
#ifndef STUDIO_BRIDGE_H
#define STUDIO_BRIDGE_H

#include "../../../../../include/marathon_spatial/marathon_spatial.h"
#include "../../../../../include/marathon_spatial/marathon_spatial_acoustics.h"
#include "../../../../../include/marathon_spatial/marathon_spatial_output.h"
#include "../../../../../include/marathon_spatial/marathon_spatial_physics.h"

#if defined(__cplusplus)
extern "C" {
#endif

/* Fills up to capacity modes of the body ms_physics_model_box builds from the same arguments, in ascending frequency.
   Returns how many were written, 0 when the library rejects the body. */
uint32_t studio_box_modes(const ms_physics_material* material, ms_vec3 size_m, float wall_thickness_m, double* frequency_hz,
                          double* damping_ratio, uint32_t capacity);

/* One room. absorption holds MS_ACOUSTIC_BANDS values per material; scattering and transmission hold one each. triangle_materials
   indexes those materials. volume, surface_area and absorption_area are the document's exact figures, so no probe has to guess them. */
typedef struct studio_room_desc {
    const float* positions;
    uint32_t vertex_count;
    const uint32_t* indices;
    uint32_t index_count;
    const uint32_t* triangle_materials;
    const float* absorption;
    const float* scattering;
    const float* transmission;
    uint32_t material_count;
    ms_vec3 bounds_min, bounds_max, probe;
    float volume, surface_area;
    const float* absorption_area;
} studio_room_desc;

/* The scene, with its room in place and its reverb prepared. metallib may be NULL, which traces on the calling thread. NULL on failure. */
ms_acoustics* studio_acoustics_create(const studio_room_desc* room, const char* metallib, float sample_rate);

/* The render side of one room's late reverb. Pass studio_reverb_bus and this as a bus processor: sends arrive on that bus, and
   leave it as the room's stereo response. It never allocates or blocks. */
typedef struct studio_reverb studio_reverb;
studio_reverb* studio_reverb_create(ms_acoustics* acoustics, uint32_t max_frames);
void studio_reverb_destroy(studio_reverb* reverb);
void studio_reverb_bus(void* reverb, float* left, float* right, uint32_t frames, uint64_t frame);

#if defined(__cplusplus)
}
#endif
#endif
