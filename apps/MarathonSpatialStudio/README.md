# Marathon Spatial Studio

A native macOS 27 app for building a room, placing sound in it and hearing the result on any headphones. It is the tool that
comes with the Marathon Spatial library. Everything it plays or exports is rendered by the library: the app itself computes no audio
other than its procedural test signals, which are fed to the library like any other clip.

Requires macOS 27 on Apple silicon (M3 or later).

## What it does

- **Room:** a plan view and a height view. Draw any simple polygon room, add blocks, and move sources and the listener.
- **Materials:** paint walls, floor, ceiling and blocks. Edit absorption band by band to enter measured data. Each material has a
  stroke pattern as well as a colour.
- **Your own audio:** drop in any file Core Audio can read. It is copied into the project, which stays self-contained. Six test signals
  are built in, so a project sounds with no imported audio.
- **Sources:** gain, loop, start time, mute and solo, directivity with a polar plot, scripted orbit or shuttle motion, and per-source Doppler,
  air absorption, near field and reflections.
- **HRTF sets:** choose from those the library ships, the analytic sphere, or your own SOFA file; each shows its licence.
- **Head tracking:** from headphones that report motion, with recentre.
- **Ambisonic beds:** any file of (order + 1)^2 ACN/SN3D channels, rotated, looped or head-locked.
- **Reflections and reverb:** the library's simulator traces early reflections, occlusion and transmission from the room you drew, for
  sources and sounding bodies alike, and runs its late reverb from the room's volume and absorption. Exports trace on the CPU, so they
  repeat. Mirror-image paths are drawn in the room by order for you to look at; they are a picture, not what is rendered. The
  inspector shows a statistical RT60 per octave band beside the library's own figure.
- **Meters:** output peak and RMS with dB readouts, plus a level per source.
- **Physics playground:** strike plates, bars and blocks with the pointer or keys 1 to 9, at any position, with chosen impulse and mallet
  hardness, or drag to scrape. See the body's modes. Record strikes onto the timeline, and export plays them.
- **3D speakers:** the binaural render, reflections and all, through the library's crosstalk canceller, for two speakers in front of
  you. Set each speaker's angle off centre and its distance; the laptop and desk presets set both. Optionally the front camera follows
  where you sit, so the effect holds as you shift; no image is kept.
- **Stereo:** the plain speaker path is an equal-power pan with distance gain. It has no HRTF and no reflections.
- **Export:** offline render to binaural or stereo WAV at 16-bit, 24-bit or 32-bit float, 44.1 to 96 kHz, reporting peak and clipped samples.
- **Presets and documents:** eight room presets, four render presets, your own of both, and `.mspatial` project packages with undo.

## Keyboard

| Key | Action |
|---|---|
| V, S, B, P, O, K | Select, add source, add block, paint, edit outline, strike (room focused) |
| Arrows | Move the selection 10 cm; with Shift, 1 m; with Option, up and down |
| Q, E | Turn the listener |
| Space or Command-R | Play or stop |
| Delete | Delete the selection |
| Command-] | Next object |
| 1 to 9 | Strike a body (playground) |
| Command-Option-1, 2 | Room, Playground |
| Command-Shift-R | Record strikes |
| Command-E | Export |

Every object in the plan is a separate VoiceOver element with move actions, and each has full numeric controls in the inspector.

## Building

    cmake -S ../.. -B ../../build -DCMAKE_BUILD_TYPE=Release && cmake --build ../../build
    Scripts/bundle.sh
    swift test

The package finds the library's static libraries through `../../build`, or in `MARATHON_SPATIAL_LIB_DIR`. Without them a plain `swift build`
still works, but playback, strikes and export stay off and the window says so. `bundle.sh` refuses that case. SwiftPM caches manifest results,
so after changing where the libraries are pass `--manifest-cache none` once.

`bundle.sh` assembles `build/Marathon Spatial Studio.app`, with the HRTF sets, the Metal ray tracing library and the licence files any copy
of the app must carry. Use that bundle: the document type is declared in its `Info.plist`, so a bare `swift run` cannot save documents.

`StudioKit` holds the document, geometry, acoustic estimates, WAV writer and the `SpatialEngine` boundary, and has no SwiftUI in it.
`SpatialEngine` is the only route to audio. `EngineFactory` chooses the engine; without the library it picks `UnavailableEngine`, which
makes no sound and says so in the window. `StudioRender` runs the library with no device, for export and for its tests; `LibraryEngine` adds
live output around it. `CStudioBridge` exposes the library's C API to Swift, and adds the modes of a box, one room's acoustic scene built
from flat arrays, and the render-thread side of its late reverb.

Clips live in a content-addressed cache in Application Support, and a project package keeps its own copies. The cache is not trimmed.
