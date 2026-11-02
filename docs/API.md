# API guide

The public C API is in [`marathon_spatial.h`](../include/marathon_spatial/marathon_spatial.h), and its comments are
the reference. Live output is in [`marathon_spatial_output.h`](../include/marathon_spatial/marathon_spatial_output.h).
The physics audio C++ API is under
[`include/marathon_spatial/physics/`](../include/marathon_spatial/physics/). This page explains how the parts fit
together.

## The model

- **Engine:** owns sources, buffers, buses, beds and one selected HRTF set. Create it with `ms_engine_desc_default`
  and `ms_engine_create`.
- **Two sides:** the *control side* is every function except render. It may lock, allocate and take time. The
  *render side* is `ms_engine_render` and any callback the engine invokes. It never locks, allocates or blocks.
- **Publishing:** `ms_*_set_*` calls change a working copy. `ms_engine_commit` publishes all of it as one finished
  snapshot at the next block boundary, so a frame's changes arrive together. Head orientation is the one exception:
  `ms_head_set_orientation` takes effect at the next block without a commit, to keep head tracking latency low.
- **Determinism:** the same calls at the same frame positions give bit-identical output, however frames are chunked
  into render calls. Live and offline rendering use the same entry point.

## Coordinates

Right handed, in metres. +X is right, +Y is up, and -Z is forward at the identity listener. Ambisonic content is
ACN/SN3D (AmbiX).

## A minimal program

```c
ms_engine_desc desc;
ms_engine_desc_default(&desc);
ms_engine* engine;
ms_engine_create(&desc, &engine);

ms_hrtf* hrtf;
ms_hrtf_create_sphere(0.0875f, desc.sample_rate, &hrtf);
ms_engine_set_hrtf(engine, hrtf);
ms_hrtf_release(hrtf);

ms_buffer buffer;
ms_buffer_load_file(engine, "sound.wav", &buffer);

ms_source source;
ms_source_create(engine, &source);
ms_source_set_buffer(engine, source, buffer, 0, 0);
ms_source_set_pose(engine, source, (ms_vec3){1, 0, -2}, (ms_vec3){0, 0, 0});
ms_source_play(engine, source, 0);
ms_engine_commit(engine);

ms_engine_render(engine, left, right, frames);
```

Every call that can fail returns an `ms_status`. `ms_status_string` names it. The [examples](../examples/README.md)
check each one.

## HRTF sets

| Function | Use |
|---|---|
| `ms_hrtf_create_sphere` | Analytic rigid sphere. Needs no data, and is a neutral fallback and a reference |
| `ms_hrtf_load` | A cooked `.mshrtf` pack |
| `ms_hrtf_import_sofa` | Your own AES69 SOFA file of type `SimpleFreeFieldHRIR`. Not for the render side |
| `ms_hrtf_save` | Writes a set as an `.mshrtf` pack |

Sets are immutable and reference counted, and one set can serve several engines. `ms_engine_set_hrtf` builds the
ambisonic decoder, which can take a second, and the change crossfades once committed.

Generic HRTFs have real limits on elevation and front-back cues. Trying several sets and using head tracking both
help.

## Sources

`ms_source_params` covers the direct path: position, velocity and Doppler, directivity, distance model, ISO 9613-1
air absorption, near-field, propagation delay, filters, three-band EQ, bus and sends. `ms_source_set_pose` updates
position alone.

| Spatial mode | Use |
|---|---|
| `MS_SPATIAL_HRTF` | Direct convolution, falling back to the ambisonic mix past `max_hrtf_voices` |
| `MS_SPATIAL_AMBISONIC` | Always the shared mix, for dense crowds and debris |
| `MS_SPATIAL_PAN` | Head-locked equal-power stereo, for interface sounds |

A source gets audio from a buffer, or from `ms_source_set_callback` on the render side. Physics audio uses the
callback. `ms_source_play` takes an absolute frame, so starts are sample accurate.

## Reflections

`ms_source_set_reflections` publishes up to `MS_MAX_REFLECTION_TAPS` early reflections for a direct voice. `delay_seconds` is
how long after the direct sound the reflection arrives, (reflected path - direct path) / speed of sound, never the absolute
path time. The same value is right with `propagation_delay` on or off. Each tap has
a delay, an arrival direction and three band gains. A path that persists must keep its slot index, because delay and
gain glide from the old value. A ray tracer publishes these, and you can also write them yourself: see
[`bring_your_own_audio.c`](../examples/bring_your_own_audio.c).

## Beds and buses

- **Beds** are higher-order ambisonic layers, decoded through the selected HRTF. A bed can be fixed to the world or
  locked to the head.
- **Buses** mix stereo into a parent, with gain and filters. `ms_bus_set_processor` adds in-place processing on the
  render side, which is where a reverb lives.
- `ms_bed_set_bus` scales a bed by a bus and its parents. Reflections and overflow voices follow their source's bus the
  same way, so muting a group mutes all of them. Sends are the host's to scale: they go to the bus you name.
- `ms_engine_set_block_callback` runs once per block on the render side, before any source is pulled, whether or not
  anything is playing. Physics audio runs from it: pass `ms_physics_block_callback` with the `ms_physics` as user.

## Physics audio

Physics audio is a separate C++ layer that needs no engine. `ModalSynth` renders bodies into any buffer, and a source
callback places the result in space.

- **Models:** `analyticModalModel` gives closed-form bars, plates and blocks. `analyzeModal` solves a closed mesh at
  cook time. Models serialize with `encodeModalModel`.
- **Excitation:** `ModalSynth::excite` strikes a body with an impulse and a contact time, which
  `hertzContactSeconds` can estimate. `ModalSynth::drive` applies a continuous force for scrapes and rolling.
- **Materials:** `materialPresets` gives plain-data starting points. A host passes its own values in the same
  structs. The figures are typical handbook values, not measurements of a sample.
- **Threading:** a `ModalSynth` has a single owner and is not thread-safe. `registerModel` allocates, so do it at load
  time. `excite`, `drive`, `stop` and `render` do not allocate.

- **Recorded layers:** `SampleLayers.hpp` plays one-shots, seamless loops and granular clouds over any loaded sample,
  with any parameter mapped onto gain, pitch, filter and grain controls, so your own recordings can follow the
  physics too.
- **Live use:** `PhysicsStream` in `Stream.hpp` takes physics facts from one control thread and hands each modal voice
  to its own spatial source. Posting and pulling never lock or allocate.

Contact and fracture helpers are in `Contact.hpp` and `Fracture.hpp`.

From C, [`marathon_spatial_physics.h`](../include/marathon_spatial/marathon_spatial_physics.h) offers the same
features through `ms_physics_create`, `ms_physics_model_box`, `ms_physics_model_mesh`, `ms_physics_impact` and
`ms_physics_contact_update`.

## Live output

Live output is a separate library target, so offline users and tests link code that cannot open a device.

```c
float rate;
ms_output_default_sample_rate(&rate);
/* Create the engine at this rate. */
ms_output* output;
ms_output_create(engine, NULL, &output);
ms_output_start(output);
```

`ms_output_create` fails with `MS_ERR_RATE_MISMATCH` if the device rate differs from the engine's.

## Telemetry

`ms_engine_get_stats` reports levels, voice counts, overs and faults. `ms_engine_get_op_counts` counts render-side
operations, which gives a budget that does not depend on the host machine's speed.
