# Examples

Every example renders offline to a WAV file. None opens an audio device, so they are safe to run anywhere and you
choose when to listen to the result.

| Example | Shows |
|---|---|
| [`orbit_offline.c`](orbit_offline.c) | The smallest complete program: one source, the built-in sphere HRTF, offline render |
| [`bring_your_own_audio.c`](bring_your_own_audio.c) | Your own sound file and optional SOFA HRTF, with early reflections written by hand |
| [`struck_plate.cpp`](struck_plate.cpp) | Physics audio alone: a closed-form plate struck with a Hertz contact time |
| [`spatial_impact.cpp`](spatial_impact.cpp) | Physics audio placed in space: a modal body as a source callback |

For sounds to try, generate the [demo content](../demo/PROVENANCE.md):

```sh
python3 demo/synthesize_demo_content.py demo/generated --verify
```
