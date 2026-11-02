/* The physics C API from plain C: everything offline, into memory. */
#include "marathon_spatial/marathon_spatial_physics.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

static int failures = 0;
#define CHECK(x) do { if (!(x)) { ++failures; printf("FAIL: %s\n", #x); } } while (0)

static double energy(const float* y, uint32_t n) {
    double e = 0;
    for (uint32_t i = 0; i < n; ++i) e += (double)y[i] * y[i];
    return e;
}

int main(void) {
    ms_physics_desc desc;
    ms_physics_desc_default(&desc);
    desc.voices = 8;
    desc.contacts = 4;
    ms_physics* physics = NULL;
    CHECK(ms_physics_create(&desc, &physics) == MS_OK);
    desc.voices = 0;
    ms_physics* none = NULL;
    CHECK(ms_physics_create(&desc, &none) == MS_ERR_INVALID_ARGUMENT);

    ms_physics_material steel, concrete;
    CHECK(ms_physics_material_preset("steel", &steel) == MS_OK);
    CHECK(ms_physics_material_preset("concrete", &concrete) == MS_OK);
    CHECK(ms_physics_material_preset("cheese", &steel) == MS_ERR_NOT_FOUND);
    CHECK(ms_physics_material_preset_count() >= 16);
    CHECK(strlen(ms_physics_material_preset_name(0)) > 0);

    ms_modal_model plate = 0;
    ms_vec3 size = {.6f, .4f, .004f};
    CHECK(ms_physics_model_box(physics, &steel, size, 0, &plate) == MS_OK);
    ms_surface surface = 0;
    CHECK(ms_physics_surface_create(physics, .001f, 3, &surface) == MS_OK);
    ms_debris_bank bank = 0;
    CHECK(ms_physics_debris_bank_create(physics, &steel, .001f, .1f, 4, &bank) == MS_OK);

    float impulse = 0, seconds = 0;
    CHECK(ms_physics_solve_impact(&steel, &concrete, 2, .5f, INFINITY, .02f, .6f, &impulse, &seconds) == MS_OK);
    CHECK(fabsf(impulse - 1.6f) < 1e-5f && seconds > 0);

    ms_physics_impact_desc hit;
    memset(&hit, 0, sizeof hit);
    hit.instance = 5; hit.model = plate; hit.impulse_ns = .1f; hit.contact_seconds = 1e-4f; hit.weight = 1;
    CHECK(ms_physics_impact(physics, &hit) == MS_OK);
    hit.instance = 0;
    CHECK(ms_physics_impact(physics, &hit) == MS_ERR_INVALID_ARGUMENT);

    ms_physics_contact_desc contact;
    ms_physics_contact_desc_default(&contact);
    contact.contact = 7; contact.model = plate; contact.surface = surface; contact.speed_m_s = 1; contact.normal_force_n = 10;
    CHECK(ms_physics_contact_update(physics, &contact) == MS_OK);

    ms_physics_fracture_desc fracture;
    memset(&fracture, 0, sizeof fracture);
    fracture.id = 9; fracture.bank = bank; fracture.parent_mass_kg = 1; fracture.energy_j = 10; fracture.height_m = 1;
    fracture.fragments = 40; fracture.seed = 2; fracture.floor = &concrete; fracture.weight = 1;
    uint32_t hits = 0;
    CHECK(ms_physics_fracture(physics, &fracture, &hits) == MS_OK && hits >= 40);

    /* One source callback per voice, as an engine pulls them. */
    float block[8][512];
    double total = 0;
    int voiced = 0;
    for (int pass = 0; pass < 100; ++pass) {
        for (uint32_t v = 0; v < 8; ++v) {
            CHECK(ms_physics_source_callback(ms_physics_voice(physics, v), block[v], 512) == 512);
            total += energy(block[v], 512);
        }
        if (pass == 1 && ms_physics_voice_instance(physics, 0) != 0) voiced = 1;
    }
    CHECK(voiced);
    CHECK(total > 0 && isfinite(total));
    /* The block clock runs from a plain ms_block_callback. */
    ms_physics_block_callback(physics, 512, 0);
    CHECK(ms_physics_source_callback(ms_physics_voice(physics, 0), block[0], 512) == 512);
    CHECK(ms_physics_contact_end(physics, 7, .01f) == MS_OK);
    CHECK(ms_physics_stop(physics, 5) == MS_OK);
    float mix[512];
    ms_physics_render_mix(physics, mix, 512);
    CHECK(isfinite(energy(mix, 512)));

    ms_sample* sample = NULL;
    CHECK(ms_sample_load("/nonexistent/file.wav", 48000, &sample) == MS_ERR_BAD_FILE);
    ms_layer_desc layer;
    ms_layer_desc_default(&layer);
    ms_layers* layers = NULL;
    CHECK(ms_layers_create(&layer, 1, 48000, 1, &layers) == MS_ERR_INVALID_ARGUMENT); /* no sample */

    ms_physics_destroy(physics);
    printf(failures ? "physics C API: %d FAILED\n" : "physics C API: all passed\n", failures);
    return failures ? 1 : 0;
}
