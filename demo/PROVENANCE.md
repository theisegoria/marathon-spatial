# Demo content provenance

Every demo sound is synthesized by `synthesize_demo_content.py`, from arithmetic only. No recording, sample
library or machine learning model is involved, and the script is dedicated to the public domain under CC0-1.0,
so anyone may use the output for any purpose. The repository stores the script, not the audio.

To regenerate the files and prove they match this table:

```sh
python3 demo/synthesize_demo_content.py demo/generated --verify
```

| File | Contents | Licence | Origin | SHA-256 |
|---|---|---|---|---|
| `click.wav` | One-sample impulse after 100 ms of silence, 1 s | CC0-1.0 | synthesized | `7676db4c5f42d1a58dc845f0321353f8d78fd77f913309967f952fa35200bb14` |
| `noise_burst.wav` | Pink noise burst, 200 ms, seed 0x5EED0001 | CC0-1.0 | synthesized | `7acf61af2ed460f66b100170dc7d84c0622bf214d85ae6ad8e3d4dcdc24c4614` |
| `struck_bar.wav` | Eight decaying modes at free-free bar ratios, 2 s | CC0-1.0 | synthesized | `f13775328120c832b9d6b9d8e63de6eda25228d41c464dfd6ad9f356e5377c6f` |
| `machine_hum.wav` | Seamless 2 s loop, harmonics of 60 Hz | CC0-1.0 | synthesized | `f0c0ee407aed105ba4db35d5122159df8a946713cec2770c1630f44f0acd9695` |
| `sweep.wav` | Exponential sine sweep, 40 Hz to 16 kHz, 3 s | CC0-1.0 | synthesized | `fba3d6736f7a109fb1688aad1e0b6dcdf4dc52ce3463a2f847d09810f97928f9` |

All files are 48 kHz, 24-bit, mono PCM, peak normalized to -6 dBFS. A file whose digest differs here must be
regenerated or explained. Add demo sounds only if they are CC0-1.0 or synthesized, and add a line to this table.
