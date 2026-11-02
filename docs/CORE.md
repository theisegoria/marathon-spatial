# The core: how the engine renders

The public contract is [`marathon_spatial.h`](../include/marathon_spatial/marathon_spatial.h). This page explains the machinery
behind it and records what has been measured.

## The boundary

The control side and the render side share no lock. Three things cross between them.

- **One snapshot** of all continuous state: listener, source parameters, reflections, bed and bus settings. `ms_engine_commit`
  copies the working state into a spare slot and publishes it with one atomic exchange. The render side takes the newest whole
  snapshot at a block boundary, so it never sees half an update and neither side ever waits.
- **One ring of commands** for structural events: play, stop, create, destroy, attach, swap HRTF. Commands are read before the
  snapshot is taken, so any command the render side sees already has its snapshot.
- **One ring of retired objects**, render to control, so nothing is freed on the render side.

Head orientation is a single packed atomic, read at each block. Telemetry is relaxed atomics. All render memory is allocated
when the engine is created.

## One block

Fixed internal blocks of 64 to 512 frames, default 128. `ms_engine_render` takes any frame count and keeps the remainder, so
output does not depend on how calls are chunked. The live output calls this same function.

1. Take commands and the snapshot; rank sources every eight blocks so the loudest hold the direct voices.
2. **Beds:** rotate into the head frame, gliding across the block, and add to the shared ambisonic mix.
3. **Each source:** read through an eight-tap Kaiser sinc at the Doppler rate, or exactly when no rate change applies, then the
   filter chain, gain ramp and sends. Then one of three paths:
   - **Direct:** partitioned convolution (below), a near-field shelf per ear and a sinc fractional delay per ear.
   - **Shared mix:** encode into the ambisonic mix. A source past `max_hrtf_voices` lands here and crosses over in one block.
   - **Pan:** head-locked equal power.
4. **Reflections:** the source's dry signal feeds three band-split delay lines, and each path encodes into the shared mix.
5. **Decode:** convolve the mix with the set's ambisonic filters into the master.
6. **Buses:** filter, gain, processor and sum, from the leaves to the master. Limiter, finite-sample guard, meters.

## Convolution

Uniformly partitioned overlap-save on `vDSP_fft_zrip`, block size B with FFTs of 2B. The direct path has no latency beyond the
interpolator's four frames, and with a 256-tap filter at B = 128 it uses two partitions.

Filters move like this. Four neighbours are found in O(1) from an octahedral grid of candidates, and weighted by modified Shepard
interpolation, with the radius taken from the fifth neighbour so a weight reaches zero exactly as its neighbour leaves the set.
The weighted minimum-phase filter is built in the time domain and transformed, only when a direction has moved by more than about
0.25 degrees. Old and new filters then crossfade for one block. The onset delays are interpolated separately and move sample by
sample, so an interaural delay glides instead of stepping.

## HRTF sets

Measured responses go through one build path, whether cooked ahead of time or imported from a user's SOFA file: each response
splits into a minimum-phase filter and an onset delay (correlation with its own minimum-phase form below 8 kHz), the bass below
150 Hz is extended, and the set is diffuse-field equalised with solid-angle weights so sets match in tone.

## Near field

Below the set's measured distance, each ear gets its own direction, which gives parallax, and a first-order shelf with gain from
a table of the rigid sphere (Duda and Martens 1998) relative to the measured distance.

## Ambisonics

ACN/SN3D up to order 7. Rotation is by Ivanic and Ruedenberg's recursion. The binaural decode is built from the selected set by
magnitude least squares: exact least squares below N c / (2 pi a), matched magnitude with free phase above it, with the phase
carried forward and advanced by the common delay so high frequencies stay aligned in time with low ones.

## Measurements

All run offline in `tests/core/CoreTests.cpp`, which links no device code. Figures are from the M5 development machine except the
counted operations, which depend on no machine. **None of this says how it sounds, and none of it is an M3 run.**

| Measure | Result |
|---|---|
| Rigid-sphere set, onset ITD against Woodworth, 0 to 90 degrees | Worst 38.6 us |
| Rigid-sphere set, ILD against the analytic sphere | Within 1.5 dB at 1, 2, 4, 8 kHz |
| Measured sets, onset ITD against Woodworth at each set's own radius | Worst 30 to 97 us |
| Near field at 21.9 cm, ILD against the analytic sphere, 200 Hz to 2 kHz | Within 0.3 dB |
| Left-right mirror error, rigid-sphere set | -94.2 dB |
| Doppler at +-34.3 m/s | 909.09 and 1111.11 Hz, as computed |
| Air absorption over 100 m against ISO 9613-1 | Worst 4.6 dB, at 16 kHz |
| Bed against direct path, order 3 | ITD within 9.5 us under 1 kHz, ILD within 4.4 dB |
| Bed against direct path, order 5 | ITD within 2.3 us under 1 kHz, ILD within 1.8 dB |
| Block sizes 64, 256, 512 against 128 | Below -129 dB |
| Render thread, allocations and blocking calls, with commits racing | 0 |
| Output hash over one second of a busy scene, four chunkings | Identical in every bit |

Counted operations per 128-frame block at 48 kHz, with a new pose for every source at each 60 Hz frame. A flop count assumes
2.5 n log2 n per real FFT, 8 per complex multiply-accumulate and 2 per real one. The M3 column assumes one performance core
sustains 15% of its 64.8 GFLOP/s peak, an assumption and not a measurement.

| Load | kflop per block | Share of one M3 performance core |
|---|---|---|
| 32 direct | 2002 | 7.7% |
| 64 direct | 4004 | 15.4% |
| 64 direct, 192 shared | 7396 | 28.5% |
| 32 direct, 12 reflections each | 4781 | 18.4% |
| 64 direct, none moving | 2327 | 9.0% |

## Known limits

- Generic HRTFs have real limits on elevation and front-back cues. Measurements prove cues are present, not that a listener
  hears them.
- `ms_engine_render` runs on one thread. Nothing divides the voices across cores yet.
- Reflection gains are three bands, and a reflection only exists on a source holding a direct voice.
- Air absorption is two shelves fitted to the ISO curve, so it drifts at very long range.
- The rigid-sphere near-field model is applied to every set, whatever head was measured.
