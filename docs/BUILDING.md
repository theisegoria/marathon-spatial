# Building

## Requirements

- macOS 27 or later on Apple M3 or later
- Xcode 27, for the macOS 27 SDK and Swift 6.2 or later
- CMake 3.28 or later

The library has no dependencies to install. Everything it needs is an Apple framework or is vendored in
`third_party/`.

## Library, tests and examples

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --label-regex cpu --output-on-failure
```

Tests render offline into memory and never open an audio device. Tests labelled `gpu` need a Metal device and run
with `ctest --test-dir build --label-regex gpu`.

The examples build into `build/examples`, and each writes a WAV file:

```sh
build/examples/orbit_offline orbit.wav
```

## Marathon Spatial Studio

Build the library first, because the app links its static libraries:

```sh
cd apps/MarathonSpatialStudio
Scripts/bundle.sh
swift test
```

`bundle.sh` stops if any library is missing. A plain `swift build` without them makes a variant that cannot sound and says so
in its window. `swift test` renders offline into memory, and its binary does not link the output library.

See the [Studio README](../apps/MarathonSpatialStudio/README.md).

## Demo content

```sh
python3 demo/synthesize_demo_content.py demo/generated --verify
```

This needs Python 3.10 or later and nothing else. The output goes to an ignored directory.

## Using it from your own project

Add this repository as a subdirectory and link the target you need:

```cmake
add_subdirectory(marathon-spatial)
target_link_libraries(your_target PRIVATE marathon_spatial_core)
```

| Target | Contains |
|---|---|
| `marathon_spatial_core` | The renderer. Cannot open an audio device |
| `marathon_spatial_output` | Live Core Audio output, plus the core |
| `marathon_spatial` | Core and output together, for a program that plays live |
| `marathon_spatial_acoustics` | Geometric acoustics and late reverb, with its Metal library `acoustic_trace.metallib` |
| `marathon_spatial_physics` | Physics audio. Needs no engine |
