#!/usr/bin/env python3
"""Rebuild data/hrtf/ from the pinned upstream archives. Explicit maintenance only; builds never download.

    python3 -I tools/vendor_hrtf.py --cache DIR --cook PATH_TO_marathon_spatial_hrtf_cook [--fetch]

Every archive is pinned by URL and SHA-256. Licences were read from the primary sources named in data/hrtf/PROVENANCE.md.
Nothing here accepts a licence that needs an account or a click.
"""
import argparse
import array
import hashlib
import math
import re
import struct
import subprocess
import urllib.request
import wave
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "data" / "hrtf"
ZENODO = "https://zenodo.org/records/12092466/files/"
SETS = [
    {"out": "sadie2_d1_ku100", "archive": "D1_HRIR_SOFA.zip", "url": ZENODO + "D1_HRIR_SOFA.zip?download=1",
     "sha256": "366321fa78f211bacc0ec6bea96701625b196a5db54ec16748b0eab0b9705f75",
     "member": "D1_HRIR_SOFA/D1_48K_24bit_256tap_FIR_SOFA.sofa", "name": "SADIE II D1 Neumann KU 100", "licence": "Apache-2.0, University of York"},
    {"out": "sadie2_d2_kemar", "archive": "D2_HRIR_SOFA.zip", "url": ZENODO + "D2_HRIR_SOFA.zip?download=1",
     "sha256": "2561af2c09688c95f3909303cd2ac51cb8c15abb574903f17f635e1376cbce94",
     "member": "D2_HRIR_SOFA/D2_48K_24bit_256tap_FIR_SOFA.sofa", "name": "SADIE II D2 KEMAR", "licence": "Apache-2.0, University of York"},
    {"out": "sadie2_h3", "archive": "H3_HRIR_SOFA.zip", "url": ZENODO + "H3_HRIR_SOFA.zip?download=1",
     "sha256": "8f0b198245e872458313492879af77d0d3ef937eaedb141a9adbe1e69c26dbda",
     "member": "H3_HRIR_SOFA/H3_48K_24bit_256tap_FIR_SOFA.sofa", "name": "SADIE II H3 human subject", "licence": "Apache-2.0, University of York"},
]
MIT = {"out": "mit_kemar", "archive": "mit_kemar_full.zip", "url": "https://sound.media.mit.edu/resources/KEMAR/full.zip",
       "sha256": "3fb533f7af95a6c07944f31a2a1d9dca0c77ab525a44865fc508009efb692ef6",
       "name": "MIT KEMAR normal pinna", "licence": "MIT Media Lab 1994, free with citation"}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch(entry, cache, allow):
    path = cache / entry["archive"]
    if not path.exists():
        if not allow:
            raise SystemExit(f"missing {path}; pass --fetch to download {entry['url']}")
        urllib.request.urlretrieve(entry["url"], path)
    digest = sha256(path)
    if digest != entry["sha256"]:
        raise SystemExit(f"{path}: sha256 {digest}, pinned {entry['sha256']}")
    return path


def mit_raw(archive, destination):
    """The left-ear measurements of the normal pinna, mirrored for the right ear: the dataset's own 'compact' convention.

    Filenames give elevation and azimuth, clockwise from the front, so the source is to the right for azimuth 0..180.
    """
    pattern = re.compile(r"elev(-?\d+)/L-?\d+e(\d+)a\.wav$")
    left = {}
    with zipfile.ZipFile(archive) as z:
        for name in z.namelist():
            m = pattern.search(name)
            if not m:
                continue
            with z.open(name) as member, wave.open(member) as w:
                if (w.getnchannels(), w.getsampwidth(), w.getframerate()) != (1, 2, 44100):
                    raise SystemExit(f"{name}: unexpected format")
                samples = array.array("h", w.readframes(w.getnframes()))
            left[(int(m.group(1)), int(m.group(2)))] = [s / 32768.0 for s in samples]
    taps = len(next(iter(left.values())))
    azimuths = {}
    for elevation, azimuth in left:
        azimuths.setdefault(elevation, []).append(azimuth)
    directions, samples = [], []
    for (elevation, azimuth), ir in sorted(left.items()):
        mirror = (360 - azimuth) % 360
        if (elevation, mirror) not in left:
            raise SystemExit(f"no mirror for {elevation} {azimuth}")
        a, e = math.radians(azimuth), math.radians(elevation)
        directions += [math.cos(e) * math.sin(a), math.sin(e), -math.cos(e) * math.cos(a)]  # Head frame.
        samples += ir + left[(elevation, mirror)]
    count = len(directions) // 3
    with open(destination, "wb") as f:
        f.write(struct.pack("<IIff", count, taps, 44100.0, 1.4))
        f.write(array.array("f", directions).tobytes())
        f.write(array.array("f", samples).tobytes())
    return count


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--cook", type=Path, required=True)
    parser.add_argument("--fetch", action="store_true")
    args = parser.parse_args()
    args.cache.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    for entry in SETS:
        archive = fetch(entry, args.cache, args.fetch)
        with zipfile.ZipFile(archive) as z:
            sofa = args.cache / Path(entry["member"]).name
            sofa.write_bytes(z.read(entry["member"]))
        target = OUT / (entry["out"] + ".mshrtf")
        subprocess.run([args.cook, "--sofa", sofa, "--out", target, "--name", entry["name"], "--licence", entry["licence"]], check=True)
        print(target.name, sha256(target))
    archive = fetch(MIT, args.cache, args.fetch)
    raw = args.cache / "mit_kemar.rawhrir"
    print("mit directions", mit_raw(archive, raw))
    target = OUT / (MIT["out"] + ".mshrtf")
    subprocess.run([args.cook, "--raw", raw, "--out", target, "--name", MIT["name"], "--licence", MIT["licence"]], check=True)
    print(target.name, sha256(target))


if __name__ == "__main__":
    main()
