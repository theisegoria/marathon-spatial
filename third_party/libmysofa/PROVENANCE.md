# libmysofa

- **Upstream:** <https://github.com/hoene/libmysofa>
- **Version:** tag `v1.3.5`, commit `6cc5b15a73e9bd97810d03767082edda7f315881`.
- **Archive:** <https://github.com/hoene/libmysofa/archive/refs/tags/v1.3.5.tar.gz>, SHA-256
  `f29508c335c83d8703f943ffc9ca783ac39aca84e851357f13a55af0f8143137`.
- **Licence:** BSD-3-Clause, read from `LICENSE` in the archive on 2026-11-02 and kept here. The exception is
  `src/hrtf/portable_endian.h`, whose header is a public domain statement by Mathias Panzenböck with a fallback to BSD, MIT or
  Apache terms. `REUSE.toml` records it as `LicenseRef-Portable-Endian`.
- **Taken:** `src/hdf/*` and, from `src/hrtf`, `reader.c`, `check.c`, `tools.c`, `spherical.c`, `tools.h`, `mysofa.h` and
  `portable_endian.h`. Upstream's tests, sample SOFA files and resampler are left out.
- **Added here:** `src/hrtf/config.h` and `src/hrtf/mysofa_export.h`, which upstream generates from CMake templates. No
  upstream file is modified.
- **Use:** only `mysofa_load`, `mysofa_tocartesian`, `mysofa_getAttribute` and `mysofa_free`, on the control side. A SOFA
  file is untrusted input: `src/hrtf/HrtfFile.cpp` checks every dimension before reading the arrays.
