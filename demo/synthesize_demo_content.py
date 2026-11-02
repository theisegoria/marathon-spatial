#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 The Marathon Spatial Authors
# SPDX-License-Identifier: CC0-1.0
"""Generate the demo sounds from arithmetic alone.

    python3 demo/synthesize_demo_content.py OUTPUT_DIRECTORY [--verify]

No recording, sample library or model goes into these files, so they carry no third party rights. They are
written to disk and never played. The repository stores this script, not its output. With --verify the digests
are compared with the table in PROVENANCE.md, which fails if the arithmetic ever changes.

Every sound is 48 kHz, 24-bit, mono PCM, peak normalized to -6 dBFS, because the renderer places mono sources.
"""
from __future__ import annotations

import argparse
import hashlib
import math
from pathlib import Path
import re
import struct
import sys
import wave

RATE = 48000
PEAK = 10 ** (-6 / 20)


class Noise:
    """A fixed generator, so output never depends on the Python version."""

    def __init__(self, seed: int) -> None:
        self.state = seed & 0xFFFFFFFF or 1

    def uniform(self) -> float:
        self.state ^= (self.state << 13) & 0xFFFFFFFF
        self.state ^= self.state >> 17
        self.state ^= (self.state << 5) & 0xFFFFFFFF
        return self.state / 0x7FFFFFFF - 1.0


def pink(noise: Noise, count: int) -> list[float]:
    """White noise through three one-pole filters approximating a 1/f spectrum."""
    b0 = b1 = b2 = 0.0
    out = []
    for _ in range(count):
        white = noise.uniform()
        b0 = 0.99765 * b0 + white * 0.0990460
        b1 = 0.96300 * b1 + white * 0.2965164
        b2 = 0.57000 * b2 + white * 1.0526913
        out.append(b0 + b1 + b2 + white * 0.1848)
    return out


def click() -> list[float]:
    """A one-sample impulse with leading silence: the clearest test of timing and reflections."""
    out = [0.0] * RATE
    out[RATE // 10] = 1.0
    return out


def noise_burst() -> list[float]:
    """200 ms of pink noise with 5 ms ramps: broadband, so localization cues have energy to act on."""
    count = RATE // 5
    ramp = RATE // 200
    source = pink(Noise(0x5EED0001), count)
    return [value * min(1.0, index / ramp, (count - 1 - index) / ramp) for index, value in enumerate(source)]


def struck_bar() -> list[float]:
    """Eight inharmonic decaying modes with free-free bar frequency ratios, from a hit at t = 0."""
    ratios = [1.0, 2.757, 5.404, 8.933, 13.344, 18.637, 24.812, 31.870]
    fundamental = 440.0
    count = RATE * 2
    out = [0.0] * count
    for index, ratio in enumerate(ratios):
        frequency = fundamental * ratio
        if frequency >= RATE / 2:
            continue
        gain = 1.0 / (index + 1)
        decay = 2.0 + 3.0 * index
        for sample in range(count):
            time = sample / RATE
            out[sample] += gain * math.exp(-decay * time) * math.sin(2 * math.pi * frequency * time)
    return out


def machine_hum() -> list[float]:
    """A two-second seamless loop: harmonics of 60 Hz, each a whole number of cycles in the loop."""
    count = RATE * 2
    return [
        sum(math.sin(2 * math.pi * 60 * harmonic * sample / RATE) / harmonic for harmonic in (1, 2, 3, 5, 7))
        for sample in range(count)
    ]


def sweep() -> list[float]:
    """Exponential sine sweep, 40 Hz to 16 kHz over three seconds, with 10 ms ramps."""
    count = RATE * 3
    ramp = RATE // 100
    low, high = 40.0, 16000.0
    rate = math.log(high / low)
    out = []
    for sample in range(count):
        time = sample / count
        phase = 2 * math.pi * low * 3.0 / rate * (math.exp(rate * time) - 1.0)
        out.append(math.sin(phase) * min(1.0, sample / ramp, (count - 1 - sample) / ramp))
    return out


SOUNDS = {
    "click.wav": click,
    "noise_burst.wav": noise_burst,
    "struck_bar.wav": struck_bar,
    "machine_hum.wav": machine_hum,
    "sweep.wav": sweep,
}


def write(path: Path, samples: list[float]) -> str:
    scale = PEAK / max(abs(value) for value in samples)
    frames = bytearray()
    for value in samples:
        frames += struct.pack("<i", round(value * scale * 8388607))[:3]
    with wave.open(str(path), "wb") as out:
        out.setnchannels(1)
        out.setsampwidth(3)
        out.setframerate(RATE)
        out.writeframes(bytes(frames))
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("output", type=Path)
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    expected = dict(re.findall(r"`([a-z_]+\.wav)`.*?`([0-9a-f]{64})`", Path(__file__).with_name("PROVENANCE.md").read_text()))
    failed = False
    for name, build in SOUNDS.items():
        digest = write(args.output / name, build())
        print(digest, name)
        if args.verify and expected.get(name) != digest:
            print(f"error: {name} does not match PROVENANCE.md", file=sys.stderr)
            failed = True
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
