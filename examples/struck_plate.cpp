// SPDX-FileCopyrightText: 2026 The Marathon Spatial Authors
// SPDX-License-Identifier: Apache-2.0
//
// Physics-based audio from a description of an object, with no recording.
//
// Builds a closed-form model of a steel plate, strikes it three times with a
// contact time estimated from Hertz theory, and writes the result to a WAV file.
// It renders offline into memory and never opens an audio device.
//
//     struck_plate OUTPUT.wav
#include "marathon_spatial/physics/Modal.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <memory>
#include <vector>

namespace physics=mspatial::physics;

namespace {
constexpr double SampleRate=48000;

// 24-bit mono PCM, the only container this example needs.
bool writeWav(const char* path,const std::vector<float>& samples) {
    FILE* file=std::fopen(path,"wb");
    if(!file) return false;
    const auto u32=[&](std::uint32_t value){std::fwrite(&value,4,1,file);};
    const auto u16=[&](std::uint16_t value){std::fwrite(&value,2,1,file);};
    const std::uint32_t dataBytes=static_cast<std::uint32_t>(samples.size()*3);
    std::fwrite("RIFF",1,4,file);u32(36+dataBytes);std::fwrite("WAVEfmt ",1,8,file);
    u32(16);u16(1);u16(1);u32(static_cast<std::uint32_t>(SampleRate));u32(static_cast<std::uint32_t>(SampleRate)*3);u16(3);u16(24);
    std::fwrite("data",1,4,file);u32(dataBytes);
    for(const float sample:samples) {
        const auto value=static_cast<std::int32_t>(std::lround(std::clamp(sample,-1.f,1.f)*8388607.f));
        const unsigned char bytes[3]{static_cast<unsigned char>(value),static_cast<unsigned char>(value>>8),static_cast<unsigned char>(value>>16)};
        std::fwrite(bytes,1,3,file);
    }
    return std::fclose(file)==0;
}
}

int main(int count,char** arguments) {
    if(count!=2) {std::fprintf(stderr,"usage: %s OUTPUT.wav\n",arguments[0]);return 2;}

    // Typical handbook figures for steel. Any host passes its own values.
    const physics::ModalElastic steel{.youngModulusPa=200e9,.densityKgPerM3=7850,.poissonRatio=.29,.lossFactor=.002};
    auto model=std::make_shared<physics::ModalModel>(physics::analyticModalModel(physics::ModalBody{
        .shape=physics::ModalShape::Plate,.size={.6f,.4f,.003f},.material=steel,.maximumModes=64}));

    physics::ModalSynth synth{physics::ModalSynthSettings{.sampleRateHz=SampleRate,.maximumModesPerVoice=64}};
    const auto handle=synth.registerModel(model);

    // A 100 g steel ball at 3 m/s on a body too massive to move: harder hits are
    // shorter, which puts more energy in the high modes.
    const double contactSeconds=physics::hertzContactSeconds(physics::HertzContact{
        .relativeNormalSpeed=3,.massAKg=.1,.massBKg=std::numeric_limits<double>::infinity(),.effectiveRadiusMeters=.015,
        .youngModulusAPa=steel.youngModulusPa,.poissonRatioA=steel.poissonRatio,
        .youngModulusBPa=steel.youngModulusPa,.poissonRatioB=steel.poissonRatio});

    constexpr std::size_t Block=1024;
    constexpr std::size_t TotalFrames=static_cast<std::size_t>(SampleRate)*3;
    const struct {std::size_t frame;physics::Vec3 position;double impulse;} hits[]{
        {0,{.1f,.05f,0},.3},{24000,{-.2f,.1f,0},.15},{48000,{0,0,0},.6}};

    std::vector<float> output(TotalFrames);
    std::size_t next=0;
    for(std::size_t frame=0;frame<TotalFrames;frame+=Block) {
        for(;next<std::size(hits)&&hits[next].frame<frame+Block;++next) synth.excite(physics::ModalExcitation{
            .instance=1,.model=handle,.localPosition=hits[next].position,.impulseNs=hits[next].impulse,.contactSeconds=contactSeconds});
        const std::size_t frames=std::min(Block,TotalFrames-frame);
        synth.render(std::span<float>{output}.subspan(frame,frames),{},frames);
    }

    // Normalize for the file. Absolute level belongs to the host's mix.
    float peak=0;
    for(const float sample:output) peak=std::max(peak,std::abs(sample));
    if(peak==0) {std::fprintf(stderr,"the plate rendered silence\n");return 1;}
    for(float& sample:output) sample*=.5f/peak;
    if(!writeWav(arguments[1],output)) {std::fprintf(stderr,"could not write %s\n",arguments[1]);return 1;}
    std::printf("contact %.1f us, %zu modes, %zu frames, wrote %s\n",contactSeconds*1e6,model->modes.size(),output.size(),arguments[1]);
    return 0;
}
