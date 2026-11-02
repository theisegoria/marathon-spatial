# HRTF sets

Each `.mshrtf` pack is a cooked derivative, not the upstream file: responses are split into minimum-phase filters and
onset delays, resampled to 48 kHz where needed, diffuse-field equalised with the bass extended below 150 Hz, and truncated
to 256 taps. `tools/vendor_hrtf.py` rebuilds them from the archives below, verifying each SHA-256 before use.

Cooked on 2026-11-02 with `tools/hrtf_cook.cpp` at the state of this directory.

| Pack | SHA-256 of pack | Directions | Source set |
|---|---|---|---|
| `sadie2_d1_ku100.mshrtf` | `2622be0a6dd6494dc409e69c0cdc04efa0e8426e955dd47fa29b6214adcc0174` | 8802 | SADIE II D1, Neumann KU 100 dummy head |
| `sadie2_d2_kemar.mshrtf` | `8103e6c34e7e1cd3367f875c39068254df48afcb0df9eaf30d8e3c482ac7228c` | 8802 | SADIE II D2, KEMAR dummy head |
| `sadie2_h3.mshrtf` | `7729537b18a76722763699abac8b5801b51a23bc6eb576ee6f2cf984db6e1d4b` | 2818 | SADIE II H3, human subject |
| `mit_kemar.mshrtf` | `5f653704320f75a586ed5113b6376266c3c7323228a1ddef90e9b84ccda702c5` | 710 | MIT KEMAR, normal pinna, left ear mirrored |

## SADIE II Database

- **Publisher:** The Audio Lab, University of York, UK.
- **Record:** <https://zenodo.org/records/12092466>, DOI 10.5281/zenodo.12092466, version 2-2. Project page:
  <https://www.york.ac.uk/sadie-project/database.html>.
- **Licence, read 2026-11-02 at both locations:** Apache-2.0. The Zenodo record lists `apache2.0`, and the project page
  states "Copyright 2018, University of York. Licensed under the Apache License, Version 2.0". No account or click needed.
- **Required notice:** the project's notice, "This product includes data developed at The Audio Lab, University of York,
  UK (Cal Armstrong, Lewis Thresh, Gavin Kearney) as part of the SADIE Project", is reproduced in the library's `NOTICE`,
  as Apache-2.0 section 4 requires. The cooked files are modified from the originals, as described above.
- **Reference:** C. Armstrong, L. Thresh, D. Murphy, G. Kearney, "A Perceptual Evaluation of Individual and Non-Individual
  HRTFs: A Case Study of the SADIE II Database", Applied Sciences 8(11), 2018.

| Archive | Member used | SHA-256 of archive | MD5 published by Zenodo |
|---|---|---|---|
| `D1_HRIR_SOFA.zip` | `D1_48K_24bit_256tap_FIR_SOFA.sofa` | `366321fa78f211bacc0ec6bea96701625b196a5db54ec16748b0eab0b9705f75` | `4850c1eb8e63e2d4f605edcdb4d5c883` |
| `D2_HRIR_SOFA.zip` | `D2_48K_24bit_256tap_FIR_SOFA.sofa` | `2561af2c09688c95f3909303cd2ac51cb8c15abb574903f17f635e1376cbce94` | `f555491d1b6b1a108bb4d26b63ad84a4` |
| `H3_HRIR_SOFA.zip` | `H3_48K_24bit_256tap_FIR_SOFA.sofa` | `8f0b198245e872458313492879af77d0d3ef937eaedb141a9adbe1e69c26dbda` | `145826ddba3db973bec5340322b27909` |

## MIT KEMAR

- **Authors:** Bill Gardner and Keith Martin, MIT Media Lab, 1994.
- **Source:** <https://sound.media.mit.edu/resources/KEMAR/full.zip>, SHA-256
  `3fb533f7af95a6c07944f31a2a1d9dca0c77ab525a44865fc508009efb692ef6`.
- **Licence, read 2026-11-02 from the archive's README at <https://sound.media.mit.edu/resources/KEMAR/README>:** "This
  data is Copyright 1994 by the MIT Media Laboratory. It is provided free with no restrictions on use, provided the
  authors are cited when the data is used in any research or commercial application." That is permissive and compatible
  with shipping beside Apache-2.0 code, provided the citation stays with the data.
- **Citation:** W. G. Gardner and K. D. Martin, "HRTF Measurements of a KEMAR Dummy-Head Microphone", MIT Media Lab
  Perceptual Computing Technical Report 280, 1994.
- **Use:** the 710 `L*` files, left ear with the normal pinna. The right ear is the same measurement at the mirrored
  azimuth, as the dataset's own compact form does.

## Not included

The rigid sphere is computed at run time and carries no third-party data.
