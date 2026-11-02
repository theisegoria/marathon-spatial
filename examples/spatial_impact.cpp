// SPDX-FileCopyrightText: 2026 The Marathon Spatial Authors
// SPDX-License-Identifier: Apache-2.0
//
// Physics audio placed in space: a modal plate is a mono source callback, so the
// sound a body makes comes from its physics and its position comes from the HRTF.
// Struck to the listener's right, it renders to a binaural WAV file offline.
// Nothing opens an audio device.
//
//     spatial_impact OUTPUT.wav
#include "marathon_spatial/marathon_spatial.h"
#include "marathon_spatial/physics/Modal.hpp"
#include "wav_writer.h"

#include <algorithm>
#include <cstdio>
#include <memory>
#include <span>
#include <vector>

namespace physics=mspatial::physics;

namespace {
constexpr std::uint32_t SampleRate=48000,Block=256,Frames=SampleRate*3;

// Runs on the render side, so it must not lock or allocate. ModalSynth::render does neither.
std::uint32_t renderPlate(void* user,float* mono,std::uint32_t frames) {
    static_cast<physics::ModalSynth*>(user)->render(std::span<float>{mono,frames},{},frames);
    return frames;
}

bool check(ms_status status,const char* what) {
    if(status!=MS_OK) std::fprintf(stderr,"%s: %s\n",what,ms_status_string(status));
    return status==MS_OK;
}
}

int main(int count,char** arguments) {
    if(count!=2) {std::fprintf(stderr,"usage: %s OUTPUT.wav\n",arguments[0]);return 2;}

    const physics::ModalElastic aluminium{.youngModulusPa=69e9,.densityKgPerM3=2700,.poissonRatio=.33,.lossFactor=.003};
    auto model=std::make_shared<physics::ModalModel>(physics::analyticModalModel(physics::ModalBody{
        .shape=physics::ModalShape::Plate,.size={.5f,.35f,.004f},.material=aluminium,.maximumModes=64}));
    physics::ModalSynth synth{physics::ModalSynthSettings{.sampleRateHz=SampleRate,.maximumModesPerVoice=64,.maximumBlockFrames=Block,.outputGain=.05}};
    const auto handle=synth.registerModel(model);

    ms_engine_desc desc;
    ms_engine_desc_default(&desc);
    desc.sample_rate=SampleRate;
    desc.block_frames=Block;
    ms_engine* engine=nullptr;
    if(!check(ms_engine_create(&desc,&engine),"ms_engine_create")) return 1;
    ms_hrtf* hrtf=nullptr;
    if(!check(ms_hrtf_create_sphere(.0875f,SampleRate,&hrtf),"ms_hrtf_create_sphere")) return 1;
    if(!check(ms_engine_set_hrtf(engine,hrtf),"ms_engine_set_hrtf")) return 1;
    ms_hrtf_release(hrtf);

    ms_source source=0;
    if(!check(ms_source_create(engine,&source),"ms_source_create")) return 1;
    ms_source_params params;
    ms_source_params_default(&params);
    params.position={2.f,0.f,-.5f};
    if(!check(ms_source_set_params(engine,source,&params),"ms_source_set_params")) return 1;
    if(!check(ms_source_set_callback(engine,source,renderPlate,&synth),"ms_source_set_callback")) return 1;
    if(!check(ms_source_play(engine,source,0),"ms_source_play")) return 1;
    if(!check(ms_engine_commit(engine),"ms_engine_commit")) return 1;

    std::vector<float> output(static_cast<std::size_t>(Frames)*2);
    for(std::uint32_t frame=0;frame<Frames;frame+=Block) {
        // A real host calls excite from a contact event. In this offline loop the
        // render happens on this thread, so ModalSynth keeps its single owner.
        if(frame%SampleRate==0) synth.excite(physics::ModalExcitation{.instance=1,.model=handle,.localPosition={.1f,.05f,0},.impulseNs=.4,.contactSeconds=3e-4});
        const std::uint32_t frames=std::min(Block,Frames-frame);
        ms_engine_render_interleaved(engine,output.data()+static_cast<std::size_t>(frame)*2,frames);
    }
    ms_engine_destroy(engine);
    if(wav_write(arguments[1],output.data(),Frames,2,SampleRate)!=0) {std::fprintf(stderr,"could not write %s\n",arguments[1]);return 1;}
    std::printf("wrote %s\n",arguments[1]);
    return 0;
}
