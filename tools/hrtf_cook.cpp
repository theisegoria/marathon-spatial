// SPDX-License-Identifier: Apache-2.0
// Cooks an HRTF set into a .mshrtf pack, from an AES69 SOFA file or from the flat form tools/vendor_hrtf.py writes.
//   hrtf_cook --sofa in.sofa --out out.mshrtf [--name N] [--licence L] [--rate 48000] [--taps 256]
//   hrtf_cook --raw in.rawhrir --out out.mshrtf --name N --licence L [--rate 48000] [--taps 256]
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#include "hrtf/HrtfSet.hpp"

using namespace mspatial;

namespace {

// count, taps (uint32), rate, distance (float), then count * 3 floats (head frame), then count * 2 * taps floats.
std::unique_ptr<hrtf::RawHrirSet> readRaw(const char* path) {
    std::FILE* f = std::fopen(path, "rb");
    if (!f) return nullptr;
    std::uint32_t head[2];
    float meta[2];
    auto raw = std::make_unique<hrtf::RawHrirSet>();
    bool ok = std::fread(head, sizeof head, 1, f) == 1 && std::fread(meta, sizeof meta, 1, f) == 1;
    if (ok) {
        raw->taps = head[1];
        raw->sampleRate = meta[0];
        raw->distance = meta[1];
        raw->directions.resize(head[0]);
        for (auto& d : raw->directions) ok = ok && std::fread(&d, 12, 1, f) == 1;
        raw->samples.resize(static_cast<std::size_t>(head[0]) * 2 * head[1]);
        ok = ok && std::fread(raw->samples.data(), sizeof(float), raw->samples.size(), f) == raw->samples.size();
    }
    std::fclose(f);
    return ok ? std::move(raw) : nullptr;
}

}  // namespace

int main(int argc, char** argv) {
    std::string sofa, rawPath, out, name, licence;
    hrtf::BuildOptions options;
    options.sampleRate = 48000.0f;
    options.taps = 256;
    for (int i = 1; i + 1 < argc; i += 2) {
        const std::string key = argv[i];
        if (key == "--sofa") sofa = argv[i + 1];
        else if (key == "--raw") rawPath = argv[i + 1];
        else if (key == "--out") out = argv[i + 1];
        else if (key == "--name") name = argv[i + 1];
        else if (key == "--licence") licence = argv[i + 1];
        else if (key == "--rate") options.sampleRate = std::stof(argv[i + 1]);
        else if (key == "--taps") options.taps = static_cast<std::uint32_t>(std::stoul(argv[i + 1]));
        else return std::fprintf(stderr, "unknown option %s\n", key.c_str()), 2;
    }
    if (out.empty() || sofa.empty() == rawPath.empty()) return std::fprintf(stderr, "need --out and one of --sofa, --raw\n"), 2;
    std::unique_ptr<hrtf::RawHrirSet> raw;
    if (!sofa.empty()) {
        ms_status status;
        raw = hrtf::readSofa(sofa.c_str(), status);
    } else {
        raw = readRaw(rawPath.c_str());
    }
    if (!raw) return std::fprintf(stderr, "cannot read input\n"), 1;
    if (!name.empty()) raw->name = name;
    if (!licence.empty()) raw->licence = licence;
    const auto set = hrtf::buildSet(*raw, options);
    if (!set) return std::fprintf(stderr, "cannot build set\n"), 1;
    const ms_status status = hrtf::saveSet(*set, out.c_str());
    std::printf("%s: %u directions, %u taps, %.0f Hz, head radius %.4f m, distance %.2f m, largest onset delay %.1f samples\n", out.c_str(),
                set->directions, set->taps, static_cast<double>(set->sampleRate), static_cast<double>(set->headRadius), static_cast<double>(set->distance),
                static_cast<double>(set->maxDelay));
    return status == MS_OK ? 0 : 1;
}
