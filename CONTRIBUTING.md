# Contributing

Thank you for helping. Marathon Spatial exists to give macOS the best spatial and physics-based audio available
on any headphones, so the most useful contributions are new features, listening reports and bug reports with a way
to reproduce them.

## Scope

- **macOS 27 or later on Apple M3 or later, only.** Use Accelerate, Core Audio and Metal directly. Please do not
  add a portability layer. Changes that serve other platforms will be declined, however good they are.
- **The library stands alone.** It includes nothing from outside this repository, apart from dependencies vendored
  inside it. Hosts pass scenes and materials in as plain data.
- **Features first, real time always.** A feature must stay real time on an M3. We do not publish performance
  comparisons with other libraries.

## Before you open a pull request

Build and run the tests as described in [BUILDING](docs/BUILDING.md), then:

- Keep the audio thread free of locks, allocation and waiting on the GPU.
- Add a test. Tests render offline into memory and check measurements.
- **Tests must never open an audio output device.** A test that plays sound is a bug, even if it passes.
- Keep one change in each pull request and explain why it is needed.
- Update the docs and the [feature status](README.md#feature-status) if behaviour changes.

Measurements prove that cues are correct. They cannot prove that something sounds good, so say which headphones
you listened on if your change affects the sound.

## Audio and data

Everything in this repository must be free to redistribute.

- **Demo audio** must be CC0-1.0 or synthesized, with a line in [`demo/PROVENANCE.md`](demo/PROVENANCE.md).
- **Do not add audio from commercial or royalty-free sound libraries.** Most of them forbid redistribution, even
  when they allow use in your own projects.
- **Vendored code and data** go in their own directory with a `PROVENANCE.md` recording the source, pinned
  revision, SHA-256 and licence, a copy of the licence, and an entry in
  [`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md). The licence must be compatible with redistribution alongside
  Apache-2.0 and must not require an account or a click-through.
- **New files** need an SPDX header or a directory listed in [`REUSE.toml`](REUSE.toml).

## Licence of contributions

Unless you state otherwise, anything you intentionally submit is licensed under the
[Apache License 2.0](LICENSE), as section 5 of that licence describes. You keep your copyright. By submitting, you
confirm that you have the right to contribute the work.

## Reporting problems

Open an issue with your Mac, macOS version, headphones, the steps that reproduce the problem and, if you can, an
offline render that shows it. For security problems, see [SECURITY](SECURITY.md).
