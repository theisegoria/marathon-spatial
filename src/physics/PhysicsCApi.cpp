#include "marathon_spatial/marathon_spatial_physics.h"
#include "marathon_spatial/physics/Fracture.hpp"
#include "marathon_spatial/physics/SampleLayers.hpp"
#include "marathon_spatial/physics/Stream.hpp"
#include <cmath>
#include <deque>
#include <limits>
#include <new>
#include <vector>

using namespace mspatial::physics;
namespace {
constexpr std::uint32_t CommandCapacity=16384;
template<class F> ms_status guarded(F&& f) noexcept {
    try{return f();}
    catch(const std::bad_alloc&){return MS_ERR_OUT_OF_MEMORY;}
    catch(const std::invalid_argument&){return MS_ERR_INVALID_ARGUMENT;}
    catch(const std::runtime_error&){return MS_ERR_BAD_FILE;}
    catch(...){return MS_ERR_INVALID_ARGUMENT;}
}
ModalElastic elastic(const ms_physics_material& m){return {m.young_modulus_pa,m.density_kg_m3,m.poisson_ratio,m.loss_factor};}
MaterialPreset preset(const ms_physics_material& m){
    MaterialPreset p;p.name="host";p.elastic=elastic(m);p.restitution=m.restitution;p.staticFriction=m.static_friction;p.dynamicFriction=m.dynamic_friction;
    p.texture={m.corner_wavelength_m,m.slide_fluctuation,m.roll_fluctuation,m.stick_stiffness_n_m};p.panelThicknessM=m.panel_thickness_m;
    return p;
}
}
struct ms_physics {
    ContactAudio audio;
    PhysicsStream stream;
    std::deque<DebrisBank> banks;
    explicit ms_physics(const ContactAudioSettings& s):audio(s),stream(audio,CommandCapacity){}
};
struct ms_sample {std::shared_ptr<const MonoSample> mono;};
struct ms_layers {LayerInstrument instrument;};

extern "C" {
void ms_physics_desc_default(ms_physics_desc* d){if(d)*d={48000,32,16,64,1024};}
ms_status ms_physics_create(const ms_physics_desc* d,ms_physics** out){
    if(!d||!out)return MS_ERR_INVALID_ARGUMENT;
    return guarded([&]{
        ContactAudioSettings s;s.synth.sampleRateHz=d->sample_rate;s.synth.maximumVoices=d->voices;s.synth.maximumModesPerVoice=d->modes_per_voice;s.synth.maximumBlockFrames=d->max_block_frames;
        s.maximumContacts=d->contacts;
        *out=new ms_physics(s);
        return MS_OK;
    });
}
void ms_physics_destroy(ms_physics* p){delete p;}
uint32_t ms_physics_material_preset_count(void){return static_cast<uint32_t>(materialPresets().size());}
const char* ms_physics_material_preset_name(uint32_t i){return i<materialPresets().size()?materialPresets()[i].name.data():nullptr;}
ms_status ms_physics_material_preset(const char* name,ms_physics_material* out){
    if(!name||!out)return MS_ERR_INVALID_ARGUMENT;
    const auto* p=findMaterialPreset(name);
    if(!p)return MS_ERR_NOT_FOUND;
    const auto& e=p->elastic;const auto& t=p->texture;
    *out={float(e.youngModulusPa),float(e.densityKgPerM3),float(e.poissonRatio),float(e.lossFactor),float(p->restitution),float(p->staticFriction),float(p->dynamicFriction),
          float(t.cornerWavelengthM),float(t.slideFluctuation),float(t.rollFluctuation),float(t.stickStiffnessNPerM),float(p->panelThicknessM)};
    return MS_OK;
}
ms_status ms_physics_model_box(ms_physics* p,const ms_physics_material* m,ms_vec3 size,float wall,ms_modal_model* out){
    if(!p||!m||!out)return MS_ERR_INVALID_ARGUMENT;
    return guarded([&]{
        *out=p->audio.registerModel(std::make_shared<const ModalModel>(analyticModalModel(modalBodyForBox(Vec3{size.x,size.y,size.z},wall,elastic(*m)))));
        return MS_OK;
    });
}
ms_status ms_physics_model_mesh(ms_physics* p,const ms_physics_material* m,const float* vertices,uint32_t vertexCount,const uint32_t* indices,uint32_t indexCount,float spacing,uint32_t modes,ms_modal_model* out){
    if(!p||!m||!vertices||!indices||!out)return MS_ERR_INVALID_ARGUMENT;
    return guarded([&]{
        std::vector<Vec3> v(vertexCount);
        for(uint32_t i=0;i<vertexCount;++i)v[i]={vertices[3*i],vertices[3*i+1],vertices[3*i+2]};
        ModalAnalysisSettings s;s.latticeSpacingMeters=spacing;s.maximumModes=modes;
        auto result=analyzeModal({v,{indices,indexCount}},modalMaterial(preset(*m)),s);
        *out=p->audio.registerModel(std::make_shared<const ModalModel>(std::move(result.model)));
        return MS_OK;
    });
}
ms_status ms_physics_model_load(ms_physics* p,const void* bytes,size_t size,ms_modal_model* out){
    if(!p||!bytes||!out)return MS_ERR_INVALID_ARGUMENT;
    return guarded([&]{
        *out=p->audio.registerModel(std::make_shared<const ModalModel>(decodeModalModel({static_cast<const unsigned char*>(bytes),size})));
        return MS_OK;
    });
}
ms_status ms_physics_surface_create(ms_physics* p,float wavelength,uint32_t seed,ms_surface* out){
    if(!p||!out)return MS_ERR_INVALID_ARGUMENT;
    return guarded([&]{*out=p->audio.registerSurface(std::make_shared<const SurfaceProfile>(wavelength,seed));return MS_OK;});
}
ms_status ms_physics_debris_bank_create(ms_physics* p,const ms_physics_material* m,float lo,float hi,uint32_t count,ms_debris_bank* out){
    if(!p||!m||!out)return MS_ERR_INVALID_ARGUMENT;
    return guarded([&]{
        p->banks.push_back(makeDebrisBank(p->audio,preset(*m),lo,hi,count));
        *out=static_cast<ms_debris_bank>(p->banks.size()-1);
        return MS_OK;
    });
}
ms_status ms_physics_solve_impact(const ms_physics_material* a,const ms_physics_material* b,float speed,float massA,float massB,float radius,float e,float* impulse,float* seconds){
    if(!a||!b||!impulse||!seconds)return MS_ERR_INVALID_ARGUMENT;
    return guarded([&]{
        const auto r=solveImpact({speed,massA,std::isinf(massB)?std::numeric_limits<double>::infinity():double(massB),radius,e,elastic(*a),elastic(*b)});
        *impulse=float(r.impulseNs);*seconds=float(r.contactSeconds);
        return MS_OK;
    });
}
ms_status ms_physics_impact(ms_physics* p,const ms_physics_impact_desc* d){
    if(!p||!d)return MS_ERR_INVALID_ARGUMENT;
    ModalExcitation e;e.instance=d->instance;e.model=d->model;e.localPosition={d->local_position.x,d->local_position.y,d->local_position.z};
    e.impulseNs=d->impulse_ns;e.contactSeconds=d->contact_seconds;e.perceptualWeight=d->weight;
    if(!p->audio.valid(e)||!std::isfinite(d->delay_seconds)||d->delay_seconds<0)return MS_ERR_INVALID_ARGUMENT;
    return p->stream.postImpact(e,d->delay_seconds)?MS_OK:MS_ERR_LIMIT_REACHED;
}
void ms_physics_contact_desc_default(ms_physics_contact_desc* d){
    if(!d)return;
    *d={};
    d->kind=MS_CONTACT_SLIDE;d->dynamic_friction=.4f;d->static_friction=.5f;d->fluctuation=.25f;d->radius_m=.02f;d->effective_modulus_pa=1e10f;d->flats=1;d->force_scale_n=1;d->ramp_seconds=.02f;d->weight=1;
}
ms_status ms_physics_contact_update(ms_physics* p,const ms_physics_contact_desc* d){
    if(!p||!d)return MS_ERR_INVALID_ARGUMENT;
    ContactUpdate u;u.contact=d->contact;u.model=d->model;u.localPosition={d->local_position.x,d->local_position.y,d->local_position.z};
    u.surface=d->surface;u.rampSeconds=d->ramp_seconds;u.perceptualWeight=d->weight;
    auto& f=u.force;f.kind=d->kind==MS_CONTACT_ROLL?ContactKind::Roll:d->kind==MS_CONTACT_SOURCE?ContactKind::Source:ContactKind::Slide;
    f.speedMPerS=d->speed_m_s;f.normalForceN=d->normal_force_n;f.dynamicFriction=d->dynamic_friction;f.staticFriction=d->static_friction;f.fluctuation=d->fluctuation;
    f.stickStiffnessNPerM=d->stick_stiffness_n_m;f.radiusM=d->radius_m;f.effectiveModulusPa=d->effective_modulus_pa;f.outOfRound=d->out_of_round;f.flats=d->flats;
    f.source=d->source?&d->source->instrument:nullptr;f.forceScaleN=d->force_scale_n;
    if(!p->audio.valid(u))return MS_ERR_INVALID_ARGUMENT;
    return p->stream.post(u)?MS_OK:MS_ERR_LIMIT_REACHED;
}
ms_status ms_physics_contact_end(ms_physics* p,uint64_t contact,float ramp){return !p?MS_ERR_INVALID_ARGUMENT:p->stream.postEnd(contact,ramp)?MS_OK:MS_ERR_LIMIT_REACHED;}
ms_status ms_physics_stop(ms_physics* p,uint64_t instance){return !p?MS_ERR_INVALID_ARGUMENT:p->stream.postStop(instance)?MS_OK:MS_ERR_LIMIT_REACHED;}
ms_status ms_physics_fracture(ms_physics* p,const ms_physics_fracture_desc* d,uint32_t* hits){
    if(!p||!d||!d->floor||d->bank>=p->banks.size())return MS_ERR_INVALID_ARGUMENT;
    return guarded([&]{
        FractureEvent e;e.id=d->id;e.parentMassKg=d->parent_mass_kg;e.energyJ=d->energy_j;e.heightM=d->height_m;e.fragments=d->fragments;e.seed=d->seed;e.floor=elastic(*d->floor);
        const auto& bank=p->banks[d->bank];
        const auto plan=planFracture(e,bank);
        uint32_t accepted=0;
        for(const auto& h:plan.hits)accepted+=p->stream.postImpact(debrisExcitation(e,h,bank,d->weight),double(d->delay_seconds)+h.timeSeconds)?1:0;
        if(hits)*hits=accepted;
        return accepted==plan.hits.size()?MS_OK:MS_ERR_LIMIT_REACHED;
    });
}
void* ms_physics_voice(ms_physics* p,uint32_t slot){
    if(!p||slot>=p->audio.synth().settings().maximumVoices)return nullptr;
    return static_cast<BlockSource*>(&p->stream.voice(slot));
}
uint32_t ms_physics_source_callback(void* user,float* mono,uint32_t frames){return blockSourceCallback(user,mono,frames);}
void ms_physics_block_callback(void* user,uint32_t frames,uint64_t frame){if(user)blockCallback(&static_cast<ms_physics*>(user)->stream,frames,frame);}
uint64_t ms_physics_voice_instance(const ms_physics* p,uint32_t slot){return p?p->stream.instanceIn(slot):0;}
void ms_physics_render_mix(ms_physics* p,float* mono,uint32_t frames){if(p&&mono)p->stream.renderBus({mono,frames});}

ms_status ms_sample_load(const char* path,float rate,ms_sample** out){
    if(!path||!out)return MS_ERR_INVALID_ARGUMENT;
    return guarded([&]{
        *out=new ms_sample{std::make_shared<const MonoSample>(toMono(loadAudioFile(path,{.targetSampleRateHz=rate})))};
        return MS_OK;
    });
}
void ms_sample_release(ms_sample* s){delete s;}
uint64_t ms_sample_frames(const ms_sample* s){return s?s->mono->frames.size():0;}
int ms_sample_has_loop(const ms_sample* s){return s&&s->mono->loop?1:0;}
void ms_layer_desc_default(ms_layer_desc* d){
    if(!d)return;
    *d={};
    const SampleLayer l;
    d->kind=MS_LAYER_ONE_SHOT;d->gain=float(l.gain);d->pitch=float(l.pitch);d->crossfade_seconds=float(l.crossfadeSeconds);d->grain_density=float(l.grainDensityPerSecond);
    d->grain_seconds=float(l.grainSeconds);d->grain_position=float(l.grainPosition);d->grain_position_jitter=float(l.grainPositionJitter);d->grain_pitch_jitter=float(l.grainPitchJitter);
}
ms_status ms_layers_create(const ms_layer_desc* d,uint32_t count,float rate,uint32_t seed,ms_layers** out){
    if(!d||!count||!out)return MS_ERR_INVALID_ARGUMENT;
    return guarded([&]{
        std::vector<SampleLayer> layers(count);
        for(uint32_t i=0;i<count;++i){
            const auto& in=d[i];auto& l=layers[i];
            if(!in.sample||in.mapping_count>MS_MAX_LAYER_MAPPINGS)return MS_ERR_INVALID_ARGUMENT;
            l.kind=in.kind==MS_LAYER_LOOP?LayerKind::Loop:in.kind==MS_LAYER_GRANULAR?LayerKind::Granular:LayerKind::OneShot;
            l.sample=in.sample->mono;l.gain=in.gain;l.pitch=in.pitch;l.lowPassHz=in.low_pass_hz;l.crossfadeSeconds=in.crossfade_seconds;
            l.grainDensityPerSecond=in.grain_density;l.grainSeconds=in.grain_seconds;l.grainPosition=in.grain_position;l.grainPositionJitter=in.grain_position_jitter;l.grainPitchJitter=in.grain_pitch_jitter;
            for(uint32_t k=0;k<in.mapping_count;++k){
                const auto& m=in.mappings[k];
                if(m.target<MS_TARGET_GAIN||m.target>MS_TARGET_GRAIN_PITCH_JITTER)return MS_ERR_INVALID_ARGUMENT;
                l.mappings.push_back({m.parameter,static_cast<LayerTarget>(m.target),{m.input_min,m.input_max,m.output_min,m.output_max,m.exponent,m.logarithmic!=0},m.smoothing_seconds});
            }
        }
        *out=new ms_layers{LayerInstrument(std::move(layers),{.sampleRateHz=rate,.seed=seed})};
        return MS_OK;
    });
}
void ms_layers_destroy(ms_layers* l){delete l;}
void ms_layers_set_parameter(ms_layers* l,uint32_t parameter,float value){if(l)l->instrument.setParameter(parameter,value);}
void ms_layers_trigger(ms_layers* l){if(l)l->instrument.trigger();}
void ms_layers_release(ms_layers* l,float seconds){if(l)l->instrument.release(seconds);}
void* ms_layers_stream(ms_layers* l){return l?static_cast<BlockSource*>(&l->instrument):nullptr;}
}
