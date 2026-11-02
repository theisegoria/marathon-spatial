# Physics audio

Sound that follows what bodies do: strikes, scrapes, rolls, breaks. You describe objects and contacts as plain data, and
the library renders mono streams for spatial sources. It also plays and shapes your own recordings. Nothing in this code
opens an audio device.

C++ headers are in `include/marathon_spatial/physics/`, in namespace `mspatial::physics`. The C API is
`include/marathon_spatial/marathon_spatial_physics.h`.

## Modal bodies (`Modal.hpp`)

A struck object rings at its modes. Each mode has a frequency, a damping and a strength at each point on the surface.

- `analyzeModal`: finite-element analysis of a closed mesh, at cook time. It takes seconds. The model encodes to bytes
  named by a SHA-256, so cooked models can be stored.
- `analyticModalModel`: closed-form free bars and simply supported plates, and an estimate for blocks. It is instant.
- `ModalSynth`: the run-time resonator bank, with voice stealing by energy and mode culling. It takes impacts and
  continuous forces.

Damping is Rayleigh, or set from a loss factor with `modalMaterialFromLossFactor`.

## Materials (`Materials.hpp`)

`materialPresets()` is a table of plain structs: elastic constants, restitution, friction and surface texture. The elastic
figures are handbook values for the bulk material. The texture figures are starting points for tuning. None of them was
measured from a sample, so copy one and change any field.

## Contacts (`Contact.hpp`)

`ContactAudio` owns a `ModalSynth` and handles contact forces.

- **Impacts:** `solveImpact` gives the impulse and the Hertz contact time. `impact` queues a hit, sample accurate, at any
  delay.
- **Sliding:** a `SurfaceProfile` is a height field fixed in space, so the same path over the same spot gives the same
  force. Force is `fluctuation x mu x N`, and the spectrum's corner follows the slip speed.
- **Stick-slip:** set `stickStiffnessNPerM` for squeal and creak. The force is a sawtooth of period
  `2 (mu_s - mu_d) N / (k v)`, fading out above the Stribeck speed.
- **Rolling:** the surface is read through the Hertz contact patch, so a heavy ball filters out more of it. `outOfRound` and
  `flats` add a beat at `flats x v / (2 pi R)`. The force can open the contact but never pulls.
- **Recordings as force:** `ContactKind::Source` drives a body's modes from any `BlockSource`, such as your own audio.

Each contact is keyed by the touch and not the body, so two touches on one body are two voices that sum linearly.

## Fracture and debris (`Fracture.hpp`)

`makeDebrisBank` builds a few fragment models. `planFracture` turns a break into a list of hits: fragment masses from a
truncated power law, speeds that hold the released kinetic energy exactly, ballistic flight, then bounces that lose speed
by the restitution. `scheduleFracture` queues them. All hits on one bank entry share a voice, so a cluster of 1800 hits
costs at most one voice per entry.

Chunks too small to ring in the audible band keep the smallest model that does, so they still tick.

## Your own audio (`AudioAsset.hpp`, `SampleLayers.hpp`)

`loadAudioFile` reads WAV, AIFF, CAF and anything else macOS decodes, with loop points from WAV `smpl`, AIFF `INST` and
`MARK`, or a CAF loop region. `toMono` averages the channels, since a spatial source is a point.

`LayerInstrument` stacks one-shots, seamless loops and granular clouds. `ParameterMapping` connects any parameter (speed,
force, a user control) to gain, pitch, low-pass, or grain density, length, position and jitter, through a curve with
smoothing.

## Live use (`Stream.hpp`)

`PhysicsStream` puts `ContactAudio` on a render thread. One control thread posts events into a lock-free queue, and each
modal voice is a `BlockSource` that feeds one spatial source:

```cpp
PhysicsStream stream(audio);
ms_source_set_callback(engine, source[slot], blockSourceCallback, &stream.voice(slot));
// Game thread, each tick
stream.post(update);
ms_source_set_pose(engine, source[slot], poseOf(stream.instanceIn(slot)), velocity);
```

Register models and surfaces before the first voice is pulled. Rendering, posting and parameter changes never lock or
allocate, and the same calls at the same frames give the same bits for any block size.

## Limits

- Finite elements are P1 tetrahedra with a dense solve, bounded at 3600 degrees of freedom. They err on the stiff side:
  the first modes of a free plate with two cells through its thickness measure 11 to 21 % above a converged Kirchhoff
  reference. Coarse lattices cost accuracy, and the analytic models are exact only for their idealizations.
- Mode gains are surface normal velocity, a proxy for radiated sound. There is no acoustic radiation solve.
- Block models are estimates.
- The presets and all perceptual balance are untested by ear.
