# Third party notices

Marathon Spatial is licensed under Apache-2.0. These components keep their own licences. Each vendored directory has
a `PROVENANCE.md` with its source, pinned revision and SHA-256, and a copy of its licence.

| Component | Location | Licence | Use |
|---|---|---|---|
| libmysofa | [`third_party/libmysofa`](third_party/libmysofa) | BSD-3-Clause | Reads AES69 SOFA files for `ms_hrtf_import_sofa` |
| portable_endian.h, Mathias Panzenböck, inside libmysofa | [`third_party/libmysofa/src/hrtf`](third_party/libmysofa/src/hrtf/portable_endian.h) | Public domain statement, [`LicenseRef-Portable-Endian`](LICENSES/LicenseRef-Portable-Endian.txt) | Byte order helpers used by libmysofa |
| SADIE II Database, University of York | [`data/hrtf`](data/hrtf/PROVENANCE.md) | Apache-2.0 | Three HRTF sets, modified as `PROVENANCE.md` describes |
| MIT KEMAR, Gardner and Martin, 1994 | [`data/hrtf`](data/hrtf/PROVENANCE.md) | MIT Media Lab terms, free with citation, [`LicenseRef-MIT-KEMAR`](LICENSES/LicenseRef-MIT-KEMAR.txt) | One HRTF set |

An app or binary built from this repository must carry `LICENSE`, `NOTICE`, this file and `LICENSES/`:
`apps/MarathonSpatialStudio/Scripts/bundle.sh` copies them into the app's resources.

## System frameworks

The library links Apple frameworks that ship with macOS, such as Accelerate, Core Audio and Metal. They are not
redistributed here.

## Demo content

Demo sounds are synthesized by a script in this repository and dedicated to the public domain under CC0-1.0. See
[`demo/PROVENANCE.md`](demo/PROVENANCE.md).
