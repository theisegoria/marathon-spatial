# Feature comparison with Steam Audio

Steam Audio is an excellent, widely shipped library and the reason this comparison is worth making. This page
compares documented features only. It makes no performance or sound quality claims in either direction: nothing
here was measured, and a table cannot say how something sounds.

Checked on 2 November 2026 against Steam Audio
[v4.8.1](https://github.com/ValveSoftware/steam-audio/releases/tag/v4.8.1) (released 11 February 2026), its public
documentation, and `phonon.h` on the `master` branch at commit `480dd64f`. "Not documented" means we did not find
the feature in those sources, which is not proof that it cannot be done. Corrections are welcome.

The right-hand column states what this project intends and whether it has landed. See the
[feature status](../README.md#feature-status) for what you can use today.

## Where Steam Audio is ahead

| Area | Steam Audio | Marathon Spatial |
|---|---|---|
| Platforms | Windows, Linux, macOS, Android, iOS ([README][sa-repo]) | macOS 27 or later on Apple M3 or later, only |
| Engines and middleware | Unity, Unreal Engine, FMOD Studio, Wwise ([README][sa-repo]) | None. A C API and Swift, so you write the integration |
| Precomputed propagation | Baked reflections and pathing on probes ([guide][sa-guide]) | Not planned for the first release |
| Maturity | Public since 2017, with a long release history ([releases][sa-releases]) | Pre-release |
| Speaker layouts | Mono to 7.1 and custom layouts ([`phonon.h`][sa-header]) | Stereo. Headphones are the focus |

## Where the two share an approach

| Feature | Steam Audio | Marathon Spatial |
|---|---|---|
| Licence | Apache-2.0 ([licence][sa-licence]) | Apache-2.0 |
| HRTF rendering | Built-in HRTF or a custom SOFA file ([HRTF][sa-hrtf]) | A choice of HRTF sets |
| Direct path | Distance attenuation, air absorption, directivity, occlusion, transmission ([simulation][sa-sim]) | The same set |
| Ambisonics | Encode, rotate, decode to binaural or speakers, with order set at creation ([decode effect][sa-ambi]) | Higher-order beds decoded through the same HRTFs |
| Reverb driven by geometry | Convolution, parametric, hybrid and TrueAudio Next ([`phonon.h`][sa-header]) | Exact early reflections, convolution from the listener's traced energy field in first-order ambisonics, and a parametric network per room carrying the tail: hybrid |
| Head tracking | The application passes listener orientation ([decode effect][sa-ambi]) | The same, and where the headphones have motion sensors the app reads them |
| Custom ray tracing | Callback scene type ([scene][sa-scene]) | Scenes passed as plain data; ray tracing behind its own interface |
| Occlusion | Raycast or volumetric ([simulation][sa-sim]) | Traced against live geometry |
| Sound around obstacles | Pathing: probe-to-probe paths from baked data ([guide][sa-guide]) | Edge diffraction found from the live meshes every tick, one or two edges, with Kurze–Anderson wedge loss per band |

## Where Marathon Spatial adds features

| Feature | Steam Audio | Marathon Spatial |
|---|---|---|
| Physics-based sound *synthesis* | Not documented. "Physics-based" there means propagation, and the [guide][sa-guide] does not cover modal synthesis, impact or contact sounds | Modal analysis of meshes, closed-form bars, plates and blocks, modal synthesis, contact-driven impacts, scrapes and rolling, fracture, material acoustics |
| Standalone native app | Not documented. Integrations live inside game editors and middleware | Marathon Spatial Studio, a native macOS 27 app |
| Bring your own audio | Your application supplies buffers | Import your own audio into the app, place it, and render a binaural result, with no engine or code |
| Hardware ray tracing on Apple GPUs | Not documented. Scene types are built-in, Embree, Radeon Rays and custom ([scene][sa-scene]); Embree is documented for x86 and x86_64, and Radeon Rays is [Windows x64 only][sa-build] | Metal 4 ray tracing for reflections |
| Reflections from deforming geometry | Instanced meshes can move, rotate and scale but [cannot deform][sa-scene] | Scenes are plain data the host can change, so doors and destruction can alter reflections |
| Doppler | Not found in `phonon.h` | Part of the direct path |
| Near-field HRTF handling | Not documented | Part of the direct path |
| Portal-coupled late reverb | Not documented for reverb. [Pathing][sa-guide] finds probe-to-probe paths around occluders from baked data | Late reverb per room, coupled through portals |
| Apple silicon as the design target | Building from source lists macOS as x64 and arm64 for Android and iOS ([build docs][sa-build]). The [4.6.0 release](https://github.com/ValveSoftware/steam-audio/releases/tag/v4.6.0) lists Apple silicon for the Wwise plugin | Apple silicon only: Accelerate, Core Audio and Metal directly, with no portability layer |

## Details worth knowing

- **Material bands.** Steam Audio materials use three bands centred on 400 Hz, 2.5 kHz and 15 kHz
  ([scene][sa-scene]). `phonon.h` also has a compile-time `IPL_ENABLE_OCTAVE_BANDS` option for 11 bands
  ([header][sa-header]).
- **SOFA support.** Steam Audio accepts a restricted subset of SOFA, including the `SimpleFreeFieldHRIR`
  convention and `Data.Delay` of 0 ([guide][sa-guide-sofa]).
- **HRTF interpolation.** Steam Audio offers nearest and bilinear ([guide][sa-guide-sofa]).
- **GPU convolution.** Steam Audio's TrueAudio Next path is [Windows x64 only][sa-build].
- **Where convolution runs.** Steam Audio can convolve per source. Here one listener-centred response serves every source in the listener's
  room, which costs one convolution however many sources there are, while each source keeps its own exact early reflections.

[sa-repo]: https://github.com/ValveSoftware/steam-audio
[sa-licence]: https://github.com/ValveSoftware/steam-audio/blob/master/LICENSE.md
[sa-releases]: https://github.com/ValveSoftware/steam-audio/releases
[sa-header]: https://github.com/ValveSoftware/steam-audio/blob/master/core/src/core/phonon.h
[sa-guide]: https://valvesoftware.github.io/steam-audio/doc/capi/guide.html
[sa-guide-sofa]: https://github.com/ValveSoftware/steam-audio/blob/master/core/doc/guide.rst
[sa-hrtf]: https://valvesoftware.github.io/steam-audio/doc/capi/hrtf.html
[sa-sim]: https://valvesoftware.github.io/steam-audio/doc/capi/simulation.html
[sa-scene]: https://valvesoftware.github.io/steam-audio/doc/capi/scene.html
[sa-ambi]: https://valvesoftware.github.io/steam-audio/doc/capi/ambisonics-decode-effect.html
[sa-build]: https://valvesoftware.github.io/steam-audio/doc/capi/build-instructions.html

Steam Audio is a trademark of Valve Corporation. This project is not affiliated with or endorsed by Valve.
