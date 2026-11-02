#pragma once
// Bring your own audio: WAV, AIFF, CAF and anything else macOS decodes, read to
// floats through ExtAudioFile (decoding only, no device is involved). Loop points
// come from the WAV smpl chunk, the AIFF INST and MARK chunks, or a CAF loop region.
#include <cstdint>
#include <filesystem>
#include <optional>
#include <span>
#include <string>
#include <vector>

namespace mspatial::physics {
struct LoopRegion {std::uint64_t startFrame{},endFrame{};}; // end is exclusive
struct AudioAsset {
    double sampleRateHz{};
    std::uint32_t channels{};
    std::vector<float> samples;        // interleaved
    std::optional<LoopRegion> loop;    // in this asset's frames
    std::string container;             // "wav", "aiff", "caf" or a four-character file type
    [[nodiscard]] std::uint64_t frames() const noexcept {return channels?samples.size()/channels:0;}
};
struct AudioLoadOptions {
    double targetSampleRateHz{};       // 0 keeps the file's rate
    std::uint64_t maximumFrames{std::uint64_t(1)<<31};
};
// Throws std::runtime_error with the path and the OS status on failure.
[[nodiscard]] AudioAsset loadAudioFile(const std::filesystem::path& path,const AudioLoadOptions& options={});
// Loop points alone, from the file's chunks. Frames are at the file's own rate.
[[nodiscard]] std::optional<LoopRegion> readLoopRegion(const std::filesystem::path& path);

struct MonoSample {
    double sampleRateHz{};
    std::vector<float> frames;
    std::optional<LoopRegion> loop;
};
// Average of the channels: a spatial source is a point.
[[nodiscard]] MonoSample toMono(const AudioAsset& asset);
}
