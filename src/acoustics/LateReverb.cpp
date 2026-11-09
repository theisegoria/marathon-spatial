#include "acoustics/LateReverb.hpp"

#include <algorithm>
#include <complex>
#include <numbers>

namespace mspatial::acoustics {
namespace {

constexpr float kPi = std::numbers::pi_v<float>;
constexpr float kDbPerNeperOfEnergy = 4.342944819f;  // 10 log10(e).
constexpr float kInvSqrtLines = 0.25f;               // 1 / sqrt(16).
constexpr float kPrototypeDb = 6.f;
constexpr float kSectionQ = 1.2f;
constexpr std::uint32_t kSleepChunks = 512;  // Longer than the longest line, so "quiet" means the whole network.
constexpr float kSleepLevel = 1e-6f;         // -120 dB.

// Section 0 is a low shelf and section 7 a high shelf, so the bands at the edges hold their decay all the way to DC and Nyquist.
constexpr float sectionFrequency(std::size_t section) {
    if (section == 0) return 88.f;
    if (section == kReverbSections - 1) return 5660.f;
    return kBandCenterHz[section];
}

struct Coefficients {
    float b0, b1, b2, a1, a2;
};

Coefficients sectionCoefficients(std::size_t section, float gainDb, float sampleRate) {
    // Audio EQ cookbook forms. cos and sin only depend on the section, but this runs rarely enough to leave them here.
    const float A = std::pow(10.f, gainDb / 40.f);
    const float w = 2.f * kPi * sectionFrequency(section) / sampleRate;
    const float cosw = std::cos(w), sinw = std::sin(w);
    float b0, b1, b2, a0, a1, a2;
    if (section == 0 || section == kReverbSections - 1) {
        const float alpha = sinw / 2.f * std::sqrt(2.f);  // Shelf slope 1.
        const float root = 2.f * std::sqrt(A) * alpha;
        if (section == 0) {
            b0 = A * ((A + 1) - (A - 1) * cosw + root);
            b1 = 2 * A * ((A - 1) - (A + 1) * cosw);
            b2 = A * ((A + 1) - (A - 1) * cosw - root);
            a0 = (A + 1) + (A - 1) * cosw + root;
            a1 = -2 * ((A - 1) + (A + 1) * cosw);
            a2 = (A + 1) + (A - 1) * cosw - root;
        } else {
            b0 = A * ((A + 1) + (A - 1) * cosw + root);
            b1 = -2 * A * ((A - 1) + (A + 1) * cosw);
            b2 = A * ((A + 1) + (A - 1) * cosw - root);
            a0 = (A + 1) - (A - 1) * cosw + root;
            a1 = 2 * ((A - 1) - (A + 1) * cosw);
            a2 = (A + 1) - (A - 1) * cosw - root;
        }
    } else {
        const float alpha = sinw / (2.f * kSectionQ);
        b0 = 1 + alpha * A;
        b1 = -2 * cosw;
        b2 = 1 - alpha * A;
        a0 = 1 + alpha / A;
        a1 = -2 * cosw;
        a2 = 1 - alpha / A;
    }
    return {b0 / a0, b1 / a0, b2 / a0, a1 / a0, a2 / a0};
}

float responseDb(const Coefficients& c, float frequency, float sampleRate) {
    const std::complex<float> z1 = std::polar(1.f, -2.f * kPi * frequency / sampleRate), z2 = z1 * z1;
    const std::complex<float> value = (c.b0 + c.b1 * z1 + c.b2 * z2) / (1.f + c.a1 * z1 + c.a2 * z2);
    return 20.f * std::log10(std::abs(value));
}

template <class Filter>
inline float runFilter(Filter& filter, float x) {
    for (std::uint32_t k = 0; k < kReverbSections; ++k) {
        const auto& s = filter.sections[k];
        const float y = s.b0 * x + filter.z1[k];
        filter.z1[k] = s.b1 * x - s.a1 * y + filter.z2[k];
        filter.z2[k] = s.b2 * x - s.a2 * y;
        x = y;
    }
    return x;
}

bool isPrime(std::uint32_t n) {
    if (n < 2) return false;
    for (std::uint32_t d = 2; d * d <= n; ++d)
        if (n % d == 0) return false;
    return true;
}

}  // namespace

void LateReverbNetwork::prepare(const LateReverbConfig& config) {
    config_ = config;
    rooms_.assign(std::min(config.rooms, kMaxRooms), Room{});
    links_.assign(std::min(config.links, kMaxRoomLinks), Link{});
    roomSend_.assign(rooms_.size(), {});
    roomSendOther_.assign(rooms_.size(), {});
    for (std::uint32_t l = 0; l < links_.size(); ++l) links_[l].phaseStep = 2.f * kPi * (3.1f + 0.37f * float(l)) / config.sampleRate;
    linkInput_.assign(rooms_.size(), {});
    splits_.assign(rooms_.size(), {});
    for (std::uint32_t k = 0; k + 1 < kCoherenceBands; ++k) {  // Butterworth sections (Q = 1/sqrt 2); two in a row make Linkwitz-Riley.
        const float w = 2.f * kPi * kCoherenceCrossovers[k] / config.sampleRate, cw = std::cos(w), alpha = std::sin(w) / (2.f * 0.70710678f);
        const float a0 = 1.f + alpha;
        splitLow_[k] = {(1.f - cw) / 2.f / a0, (1.f - cw) / a0, (1.f - cw) / 2.f / a0, -2.f * cw / a0, (1.f - alpha) / a0};
        splitHigh_[k] = {(1.f + cw) / 2.f / a0, -(1.f + cw) / a0, (1.f + cw) / 2.f / a0, -2.f * cw / a0, (1.f - alpha) / a0};
    }

    std::size_t memory = 0;
    for (std::uint32_t r = 0; r < rooms_.size(); ++r) {
        const float mean = r < config.meanDelaySeconds.size() && config.meanDelaySeconds[r] > 0.f ? config.meanDelaySeconds[r] : 0.03f;
        const float meanSamples = std::clamp(mean, 0.008f, 0.12f) * config.sampleRate;
        Room& room = rooms_[r];
        // Prime lengths spread over a 1:3 geometric range, nudged per room so two rooms never ring alike.
        for (std::uint32_t l = 0; l < kReverbLines; ++l) {
            const float position = float(l) / float(kReverbLines - 1);
            std::uint32_t length = std::uint32_t(meanSamples * 0.55f * std::pow(3.f, position)) + r * 7 + 1;
            while (!isPrime(length)) ++length;
            room.lengths[l] = length;
            room.offsets[l] = std::uint32_t(memory);
            memory += length;
            room.totalDelay += float(length);
        }
    }
    delayMemory_.assign(memory, 0.f);
    hybrid_ = Hybrid{};
    hybridMemory_.assign(config.maxHybridDelay > 0 ? config.maxHybridDelay + 1 : 0, 0.f);

    // Invert the band interaction: the response of every section at every band center, per dB of section gain.
    double matrix[kBandCount][kBandCount * 2] = {};
    for (std::size_t s = 0; s < kReverbSections; ++s) {
        const Coefficients c = sectionCoefficients(s, kPrototypeDb, config.sampleRate);
        for (std::size_t b = 0; b < kBandCount; ++b) matrix[b][s] = responseDb(c, kBandCenterHz[b], config.sampleRate) / kPrototypeDb;
    }
    for (std::size_t b = 0; b < kBandCount; ++b) matrix[b][kBandCount + b] = 1.0;
    for (std::size_t k = 0; k < kBandCount; ++k) {
        std::size_t pivot = k;
        for (std::size_t i = k + 1; i < kBandCount; ++i)
            if (std::fabs(matrix[i][k]) > std::fabs(matrix[pivot][k])) pivot = i;
        for (std::size_t j = 0; j < kBandCount * 2; ++j) std::swap(matrix[k][j], matrix[pivot][j]);
        const double scale = 1.0 / matrix[k][k];
        for (std::size_t j = 0; j < kBandCount * 2; ++j) matrix[k][j] *= scale;
        for (std::size_t i = 0; i < kBandCount; ++i) {
            if (i == k) continue;
            const double f = matrix[i][k];
            for (std::size_t j = 0; j < kBandCount * 2; ++j) matrix[i][j] -= f * matrix[k][j];
        }
    }
    for (std::size_t i = 0; i < kBandCount; ++i)
        for (std::size_t j = 0; j < kBandCount; ++j) inverseResponse_[i][j] = float(matrix[i][kBandCount + j]);
}

void LateReverbNetwork::reset() {
    std::fill(delayMemory_.begin(), delayMemory_.end(), 0.f);
    for (auto& send : roomSend_) send.fill(0.f);
    for (auto& send : roomSendOther_) send.fill(0.f);
    for (auto& input : linkInput_) input.fill(0.f);
    for (Room& room : rooms_) {
        room.positions = {};
        room.designed = false;
        room.silentChunks = 0;
        for (Filter& filter : room.filters) {
            filter.z1 = {};
            filter.z2 = {};
        }
    }
    for (Link& link : links_) {
        link.filter.z1 = {};
        link.filter.z2 = {};
    }
    std::fill(hybridMemory_.begin(), hybridMemory_.end(), 0.f);
    hybrid_.filter.z1 = {};
    hybrid_.filter.z2 = {};
}

void LateReverbNetwork::setHybrid(std::uint32_t room, std::uint32_t delaySamples, const Bands& gainDb) {
    if (hybridMemory_.empty() || room >= rooms_.size()) room = kNone;
    delaySamples = std::min<std::uint32_t>(delaySamples, hybridMemory_.empty() ? 0u : std::uint32_t(hybridMemory_.size() - 1));
    if (room != hybrid_.room || delaySamples != hybrid_.delay) {
        std::fill(hybridMemory_.begin(), hybridMemory_.end(), 0.f);
        hybrid_.write = 0;
        hybrid_.filter.z1 = {};
        hybrid_.filter.z2 = {};
    }
    hybrid_.room = room;
    hybrid_.delay = delaySamples;
    if (room == kNone || gainDb == hybrid_.gainDb) return;
    hybrid_.gainDb = gainDb;
    float mean = 0.f, spread = 0.f;
    for (float db : gainDb) mean += db / float(kBandCount);
    for (float db : gainDb) spread = std::max(spread, std::fabs(db - mean));
    hybrid_.flatGain = std::pow(10.f, mean / 20.f);
    hybrid_.shaped = spread > 0.1f;
    if (!hybrid_.shaped) return;
    for (std::uint32_t s = 0; s < kReverbSections; ++s) {
        float gain = 0.f;
        for (std::size_t b = 0; b < kBandCount; ++b) gain += inverseResponse_[s][b] * (gainDb[b] - mean);
        const Coefficients c = sectionCoefficients(s, std::clamp(gain, -30.f, 30.f), config_.sampleRate);
        hybrid_.filter.sections[s] = {c.b0, c.b1, c.b2, c.a1, c.a2};
    }
}

std::uint32_t LateReverbNetwork::totalDelay(std::uint32_t room) const { return std::uint32_t(rooms_[room].totalDelay); }

void LateReverbNetwork::design(Room& room) {
    float mean = 0.f;
    for (float db : room.currentDbPerSample) mean += db / float(kBandCount);
    std::array<float, kBandCount> shape{};
    for (std::size_t s = 0; s < kBandCount; ++s)
        for (std::size_t b = 0; b < kBandCount; ++b) shape[s] += inverseResponse_[s][b] * (room.currentDbPerSample[b] - mean);
    for (std::uint32_t l = 0; l < kReverbLines; ++l) {
        const float length = float(room.lengths[l]);
        room.flatGain[l] = std::pow(10.f, mean * length / 20.f);
        for (std::uint32_t s = 0; s < kReverbSections; ++s) {
            const Coefficients c = sectionCoefficients(s, std::clamp(shape[s] * length, -30.f, 30.f), config_.sampleRate);
            room.filters[l].sections[s] = {c.b0, c.b1, c.b2, c.a1, c.a2};
        }
    }
}

void LateReverbNetwork::designLink(Link& link, const Bands& rate, float fromVolume, float toVolume) {
    // The link's broadband gain sits at the band mean in dB; the filter carries the deviation, which is zero for an open doorway.
    float mean = 0.f;
    std::array<float, kBandCount> db{};
    for (std::size_t b = 0; b < kBandCount; ++b) {
        db[b] = 10.f * std::log10(std::max(rate[b], 1e-9f));
        mean += db[b] / float(kBandCount);
    }
    float spread = 0.f;
    for (std::size_t b = 0; b < kBandCount; ++b) spread = std::max(spread, std::fabs(db[b] - mean));
    // Output power stands for energy density in every room. A density w in the sender raises the receiver's density by rate * w * V_from / V_to each
    // second, and the receiver stores density at totalDelay per unit of output power.
    const Room& to = rooms_[link.to];
    link.gain = std::sqrt(std::pow(10.f, mean / 10.f) * (fromVolume / toVolume) * to.totalDelay / config_.sampleRate);
    link.shaped = spread > 0.1f;
    if (!link.shaped) return;
    for (std::uint32_t s = 0; s < kReverbSections; ++s) {
        float gain = 0.f;
        for (std::size_t b = 0; b < kBandCount; ++b) gain += inverseResponse_[s][b] * (db[b] - mean);
        const Coefficients c = sectionCoefficients(s, std::clamp(gain, -30.f, 30.f), config_.sampleRate);
        link.filter.sections[s] = {c.b0, c.b1, c.b2, c.a1, c.a2};
    }
}

void LateReverbNetwork::apply(const AcousticSnapshot& snapshot) {
    for (std::uint32_t r = 0; r < rooms_.size(); ++r) {
        Room& room = rooms_[r];
        room.active = r < snapshot.roomCount && snapshot.rooms[r].active;
        if (!room.active) continue;
        const RoomAcoustics& source = snapshot.rooms[r];
        for (std::size_t b = 0; b < kBandCount; ++b) room.targetDbPerSample[b] = -kDbPerNeperOfEnergy * source.decayRate[b] / config_.sampleRate;
        // Stored energy per sample of input is g^2; the room adds Q / V of energy density, and 1 m of free field holds Q / (4 pi c).
        room.inputGain = std::sqrt(4.f * kPi * config_.speedOfSound * room.totalDelay / (source.volume * config_.sampleRate));
        if (!room.designed) {
            room.currentDbPerSample = room.targetDbPerSample;
            design(room);
            room.designed = true;
        }
    }
    for (std::uint32_t l = 0; l < links_.size(); ++l) {
        Link& link = links_[l];
        link.active = false;
        if (l >= snapshot.linkCount) continue;
        const RoomLink& source = snapshot.links[l];
        if (source.from >= rooms_.size() || source.to >= rooms_.size() || !rooms_[source.from].active || !rooms_[source.to].active) continue;
        link.active = true;
        link.from = source.from;
        link.to = source.to;
        designLink(link, source.rate, snapshot.rooms[source.from].volume, snapshot.rooms[source.to].volume);
    }
}

void LateReverbNetwork::setCoherence(const float* coherence) {
    if (!coherence) {
        binaural_.store(false, std::memory_order_release);
        return;
    }
    for (std::uint32_t k = 0; k < kCoherenceBands; ++k)
        coherence_[k].store(std::isfinite(coherence[k]) ? std::clamp(coherence[k], -1.f, 1.f) : 0.f, std::memory_order_relaxed);
    binaural_.store(true, std::memory_order_release);
}

void LateReverbNetwork::binaural(std::uint32_t index, float* left, float* right, std::uint32_t frames) {
    float a[kCoherenceBands], b[kCoherenceBands];
    for (std::uint32_t k = 0; k < kCoherenceBands; ++k) {
        const float c = coherence_[k].load(std::memory_order_relaxed);
        a[k] = std::sqrt(0.5f * (1.f + c));
        b[k] = std::sqrt(0.5f * (1.f - c));
    }
    const auto section = [](const Biquad& q, Section& s, float x) {
        const float y = q.b0 * x + s.z1;
        s.z1 = q.b1 * x - q.a1 * y + s.z2;
        s.z2 = q.b2 * x - q.a2 * y;
        return y;
    };
    std::array<Split, 2>& split = splits_[index];
    for (std::uint32_t n = 0; n < frames; ++n) {
        float shaped[2];
        const float input[2] = {left[n], right[n]};
        const float* gains[2] = {a, b};
        for (std::uint32_t c = 0; c < 2; ++c) {
            float rest = input[c], sum = 0.f;
            for (std::uint32_t k = 0; k + 1 < kCoherenceBands; ++k) {
                auto& st = split[c].stages[k];
                const float low = section(splitLow_[k], st[1], section(splitLow_[k], st[0], rest));
                rest = section(splitHigh_[k], st[3], section(splitHigh_[k], st[2], rest));
                sum += gains[c][k] * low;
            }
            shaped[c] = sum + gains[c][kCoherenceBands - 1] * rest;
        }
        left[n] = shaped[0] + shaped[1];
        right[n] = shaped[0] - shaped[1];
    }
}

void LateReverbNetwork::renderRoom(std::uint32_t index, const float* input, std::uint32_t offset, std::uint32_t frames, float* left, float* right) {
    Room& room = rooms_[index];
    float* const memory = delayMemory_.data();
    const float inputGain = room.inputGain;
    const std::array<float, kReverbChunk>& injected = linkInput_[index];
    std::array<float, kReverbChunk>& send = roomSend_[index];
    std::array<float, kReverbChunk>& sendOther = roomSendOther_[index];
    for (std::uint32_t n = 0; n < frames; ++n) {
        float x[kReverbLines];
        for (std::uint32_t l = 0; l < kReverbLines; ++l) {
            const float tap = memory[room.offsets[l] + room.positions[l]];
            x[l] = runFilter(room.filters[l], tap * room.flatGain[l]);
        }
        float sum = 0.f, patternA = 0.f, patternB = 0.f;
        for (std::uint32_t l = 0; l < kReverbLines; ++l) {
            sum += x[l];
            patternA += (l & 1) ? -x[l] : x[l];
            patternB += (l & 2) ? -x[l] : x[l];
        }
        patternA *= kInvSqrtLines;
        patternB *= kInvSqrtLines;
        if (left) left[n] = patternA;
        if (right) right[n] = patternB;
        send[n] = patternA;
        sendOther[n] = patternB;
        const float excitation = ((input ? input[offset + n] * inputGain : 0.f) + injected[n]) * kInvSqrtLines;
        const float reflect = sum * (2.f / float(kReverbLines));
        for (std::uint32_t l = 0; l < kReverbLines; ++l) {
            memory[room.offsets[l] + room.positions[l]] = x[l] - reflect + excitation;
            if (++room.positions[l] == room.lengths[l]) room.positions[l] = 0;
        }
    }
}

void LateReverbNetwork::process(const float* const* inputs, float* const* outputsLeft, float* const* outputsRight, std::uint32_t frames) {
    std::uint32_t done = 0;
    while (done < frames) {
        const std::uint32_t chunk = std::min(frames - done, kReverbChunk);
        for (std::uint32_t r = 0; r < rooms_.size(); ++r) {
            Room& room = rooms_[r];
            float* left = outputsLeft && outputsLeft[r] ? outputsLeft[r] + done : nullptr;
            float* right = outputsRight && outputsRight[r] ? outputsRight[r] + done : nullptr;
            if (!room.active) {
                if (left) std::fill_n(left, chunk, 0.f);
                if (right) std::fill_n(right, chunk, 0.f);
                continue;
            }
            bool changed = false;
            for (std::size_t b = 0; b < kBandCount; ++b) {
                const float delta = room.targetDbPerSample[b] - room.currentDbPerSample[b];
                if (std::fabs(delta) <= 2e-3f * std::fabs(room.targetDbPerSample[b])) {
                    if (delta != 0.f) {
                        room.currentDbPerSample[b] = room.targetDbPerSample[b];
                        changed = true;
                    }
                } else {
                    room.currentDbPerSample[b] += delta * 0.04f;  // About 15 ms to settle.
                    changed = true;
                }
            }
            if (changed) design(room);
            // A room with nothing in it costs nothing: once its tail has died it sleeps until something arrives.
            const float* input = inputs ? inputs[r] : nullptr;
            std::uint32_t inputOffset = done;
            if (r == hybrid_.room) {
                // The room's own send, delayed and shaped. Its delay line runs even while the network sleeps, so nothing is lost.
                const std::uint32_t size = std::uint32_t(hybridMemory_.size());
                for (std::uint32_t n = 0; n < chunk; ++n) {
                    hybridMemory_[hybrid_.write] = input ? input[done + n] : 0.f;
                    const std::uint32_t read = (hybrid_.write + size - hybrid_.delay) % size;
                    const float value = hybridMemory_[read] * hybrid_.flatGain;
                    hybrid_.chunk[n] = hybrid_.shaped ? runFilter(hybrid_.filter, value) : value;
                    hybrid_.write = (hybrid_.write + 1) % size;
                }
                input = hybrid_.chunk.data();
                inputOffset = 0;
            }
            bool driven = false;
            for (std::uint32_t n = 0; n < chunk && !driven; ++n) driven = (input && input[inputOffset + n] != 0.f) || linkInput_[r][n] != 0.f;
            if (!driven && room.silentChunks >= kSleepChunks) {
                if (left) std::fill_n(left, chunk, 0.f);
                if (right) std::fill_n(right, chunk, 0.f);
                roomSend_[r].fill(0.f);
                roomSendOther_[r].fill(0.f);
                continue;
            }
            renderRoom(r, input, inputOffset, chunk, left, right);
            if (left && right && binaural_.load(std::memory_order_acquire)) binaural(r, left, right, chunk);
            float level = 0.f;
            for (std::uint32_t n = 0; n < chunk; ++n) level = std::max(level, std::fabs(roomSend_[r][n]));
            room.silentChunks = driven || level > kSleepLevel ? 0 : room.silentChunks + 1;
        }
        for (auto& input : linkInput_) input.fill(0.f);
        for (Link& link : links_) {
            if (!link.active) continue;
            // A fixed real coupling splits the modes of the two networks: some decay slower than the diffuse field and take over the tail. Slowly
            // rotating between two uncorrelated output patterns averages that splitting away and leaves the energy exchange the theory asks for.
            link.phase += link.phaseStep * float(chunk);
            if (link.phase > 2.f * kPi) link.phase -= 2.f * kPi;
            const float c = std::cos(link.phase) * link.gain, s = std::sin(link.phase) * link.gain;
            const std::array<float, kReverbChunk>& send = roomSend_[link.from];
            const std::array<float, kReverbChunk>& sendOther = roomSendOther_[link.from];
            std::array<float, kReverbChunk>& input = linkInput_[link.to];
            for (std::uint32_t n = 0; n < chunk; ++n) {
                const float value = send[n] * c + sendOther[n] * s;
                input[n] += link.shaped ? runFilter(link.filter, value) : value;
            }
        }
        done += chunk;
    }
}

}  // namespace mspatial::acoustics
