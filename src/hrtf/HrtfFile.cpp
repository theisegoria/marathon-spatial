// SPDX-License-Identifier: Apache-2.0
// The cooked .mshrtf pack, and import of AES69 SOFA files.
#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <vector>

#include "HrtfSet.hpp"
#include "mysofa.h"

namespace mspatial::hrtf {

namespace {

struct Header {
    char magic[8];  // "MSHRTF1\n"
    std::uint32_t headerBytes, directions, taps, reserved;
    float sampleRate, headRadius, distance, maxDelay;
    char name[64], licence[64];
    std::uint64_t dataBytes, checksum;  // FNV-1a over the payload
};
constexpr char kMagic[8] = {'M', 'S', 'H', 'R', 'T', 'F', '1', '\n'};

std::uint64_t fnv(const void* data, std::size_t bytes, std::uint64_t hash) {
    const auto* p = static_cast<const unsigned char*>(data);
    for (std::size_t i = 0; i < bytes; ++i) hash = (hash ^ p[i]) * 1099511628211ull;
    return hash;
}

struct CloseFile {
    void operator()(std::FILE* f) const { std::fclose(f); }
};

}  // namespace

ms_status saveSet(const HrtfSet& set, const char* path) {
    Header h{};
    std::copy(kMagic, kMagic + 8, h.magic);
    h.headerBytes = sizeof(Header);
    h.directions = set.directions;
    h.taps = set.taps;
    h.sampleRate = set.sampleRate;
    h.headRadius = set.headRadius;
    h.distance = set.distance;
    h.maxDelay = set.maxDelay;
    std::snprintf(h.name, sizeof h.name, "%s", set.name.c_str());
    std::snprintf(h.licence, sizeof h.licence, "%s", set.licence.c_str());
    const std::size_t dirBytes = set.directions * 12ull, delayBytes = set.directions * 8ull, filterBytes = set.directions * 8ull * set.taps;
    h.dataBytes = dirBytes + delayBytes + filterBytes;
    h.checksum = fnv(set.filters.data(), filterBytes, fnv(set.delays.data(), delayBytes, fnv(set.dirs.data(), dirBytes, 14695981039346656037ull)));
    std::unique_ptr<std::FILE, CloseFile> f(std::fopen(path, "wb"));
    if (!f) return MS_ERR_NOT_FOUND;
    if (std::fwrite(&h, sizeof h, 1, f.get()) != 1 || std::fwrite(set.dirs.data(), 1, dirBytes, f.get()) != dirBytes ||
        std::fwrite(set.delays.data(), 1, delayBytes, f.get()) != delayBytes ||
        std::fwrite(set.filters.data(), 1, filterBytes, f.get()) != filterBytes)
        return MS_ERR_BAD_FILE;
    return MS_OK;
}

std::unique_ptr<HrtfSet> loadSet(const char* path, ms_status& status) {
    status = MS_ERR_NOT_FOUND;
    std::unique_ptr<std::FILE, CloseFile> f(std::fopen(path, "rb"));
    if (!f) return nullptr;
    status = MS_ERR_BAD_FILE;
    Header h{};
    if (std::fread(&h, sizeof h, 1, f.get()) != 1 || !std::equal(kMagic, kMagic + 8, h.magic) || h.headerBytes != sizeof(Header)) return nullptr;
    if (h.directions == 0 || h.directions > 200000 || h.taps == 0 || h.taps > 512 || !(h.sampleRate >= 8000.0f && h.sampleRate <= 384000.0f)) return nullptr;
    const std::size_t dirBytes = h.directions * 12ull, delayBytes = h.directions * 8ull, filterBytes = h.directions * 8ull * h.taps;
    if (h.dataBytes != dirBytes + delayBytes + filterBytes) return nullptr;
    auto set = std::make_unique<HrtfSet>();
    h.name[sizeof h.name - 1] = 0;
    h.licence[sizeof h.licence - 1] = 0;
    set->name = h.name;
    set->licence = h.licence;
    set->sampleRate = h.sampleRate;
    set->headRadius = h.headRadius;
    set->distance = h.distance;
    set->maxDelay = h.maxDelay;
    set->directions = h.directions;
    set->taps = h.taps;
    set->dirs.reset(h.directions * 3ull);
    set->delays.reset(h.directions * 2ull);
    set->filters.reset(h.directions * 2ull * h.taps);
    if (std::fread(set->dirs.data(), 1, dirBytes, f.get()) != dirBytes || std::fread(set->delays.data(), 1, delayBytes, f.get()) != delayBytes ||
        std::fread(set->filters.data(), 1, filterBytes, f.get()) != filterBytes)
        return nullptr;
    if (fnv(set->filters.data(), filterBytes, fnv(set->delays.data(), delayBytes, fnv(set->dirs.data(), dirBytes, 14695981039346656037ull))) != h.checksum) return nullptr;
    if (!(set->headRadius > 0.03f && set->headRadius < 0.2f) || !(set->distance > 0.0f) || !(set->maxDelay >= 0.0f)) return nullptr;
    set->finalize();
    status = MS_OK;
    return set;
}

bool receiversSwapped(const float* responses, const float* positions, std::uint32_t measurements, std::uint32_t taps, bool metadataSwapped) {
    double evidence = 0.0, total = 0.0;  // Positive when the first receiver is louder for sources on the +y side, the left.
    for (std::uint32_t m = 0; m < measurements; ++m) {
        const float* p = positions + m * 3ull;
        const double radius = std::sqrt(static_cast<double>(p[0]) * p[0] + static_cast<double>(p[1]) * p[1] + static_cast<double>(p[2]) * p[2]);
        if (!(radius > 0.0) || std::fabs(p[1]) < 0.5 * radius) continue;
        const float* first = responses + m * 2ull * taps;
        const float* second = first + taps;
        double energy[2] = {0.0, 0.0};
        for (std::uint32_t n = 0; n < taps; ++n) {
            energy[0] += static_cast<double>(first[n]) * first[n];
            energy[1] += static_cast<double>(second[n]) * second[n];
        }
        evidence += (energy[0] - energy[1]) * (p[1] > 0.0f ? 1.0 : -1.0);
        total += energy[0] + energy[1];
    }
    return total > 0.0 && std::fabs(evidence) > 0.05 * total ? evidence < 0.0 : metadataSwapped;
}

std::unique_ptr<RawHrirSet> readSofa(const char* path, ms_status& status) {
    status = MS_ERR_BAD_FILE;
    int error = 0;
    MYSOFA_HRTF* sofa = mysofa_load(path, &error);
    if (!sofa) {
        status = error == MYSOFA_OK ? MS_ERR_BAD_FILE : (error == ENOENT || error == EIO ? MS_ERR_NOT_FOUND : MS_ERR_BAD_FILE);
        return nullptr;
    }
    struct Free {
        void operator()(MYSOFA_HRTF* s) const { mysofa_free(s); }
    };
    std::unique_ptr<MYSOFA_HRTF, Free> guard(sofa);
    const std::uint32_t M = sofa->M, R = sofa->R, N = sofa->N;
    if (M == 0 || M > 200000 || R != 2 || N < 8 || N > 65536) return nullptr;
    if (sofa->DataIR.elements != M * R * N || !sofa->DataIR.values || !sofa->SourcePosition.values || sofa->SourcePosition.elements < M * 3 ||
        !sofa->DataSamplingRate.values || sofa->DataSamplingRate.elements < 1)
        return nullptr;
    const float rate = sofa->DataSamplingRate.values[0];
    if (!(rate >= 8000.0f && rate <= 384000.0f)) return nullptr;
    mysofa_tocartesian(sofa);  // SOFA: +x front, +y left, +z up. Source positions relative to the listener.

    bool metadataSwapped = false;
    if (sofa->ReceiverPosition.values && sofa->ReceiverPosition.elements >= 6)
        metadataSwapped = sofa->ReceiverPosition.values[1] < sofa->ReceiverPosition.values[sofa->ReceiverPosition.elements / 2 + 1];
    const bool swapEars = receiversSwapped(sofa->DataIR.values, sofa->SourcePosition.values, M, N, metadataSwapped);

    auto raw = std::make_unique<RawHrirSet>();
    raw->taps = N;
    raw->sampleRate = rate;
    raw->directions.resize(M);
    raw->samples.resize(static_cast<std::size_t>(M) * 2 * N);
    double distance = 0.0;
    for (std::uint32_t m = 0; m < M; ++m) {
        const float* p = sofa->SourcePosition.values + m * 3ull;
        const Vec3 position{-p[1], p[2], -p[0]};
        distance += core::length(position);
        raw->directions[m] = core::normalized(position);
        for (std::uint32_t ear = 0; ear < 2; ++ear) {
            const float* in = sofa->DataIR.values + (static_cast<std::size_t>(m) * 2 + (swapEars ? 1 - ear : ear)) * N;
            float* out = raw->samples.data() + (static_cast<std::size_t>(m) * 2 + ear) * N;
            std::copy(in, in + N, out);
        }
    }
    raw->distance = static_cast<float>(distance / M);
    // Data.Delay holds onsets the responses themselves leave out, which in some files is the entire ITD. AES69 allows one pair
    // for the whole set or one pair per measurement.
    const std::uint32_t storedDelays = sofa->DataDelay.values ? sofa->DataDelay.elements : 0;
    if (storedDelays == 2 || storedDelays == M * 2) {
        raw->delays.resize(static_cast<std::size_t>(M) * 2);
        const float* stored = sofa->DataDelay.values;
        for (std::uint32_t m = 0; m < M; ++m)
            for (std::uint32_t ear = 0; ear < 2; ++ear) {
                const float value = stored[(storedDelays == 2 ? 0 : m * 2ull) + (swapEars ? 1 - ear : ear)];
                if (!std::isfinite(value) || value < 0.0f) return nullptr;
                raw->delays[m * 2ull + ear] = value;
            }
    }
    if (char* title = mysofa_getAttribute(sofa->attributes, const_cast<char*>("Title"))) raw->name = title;
    if (char* licence = mysofa_getAttribute(sofa->attributes, const_cast<char*>("License"))) raw->licence = licence;
    status = MS_OK;
    return raw;
}

}  // namespace mspatial::hrtf
