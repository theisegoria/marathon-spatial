#include "marathon_spatial/physics/SampleLayers.hpp"
#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <string>

namespace mspatial::physics {
namespace {
constexpr double Pi=std::numbers::pi;
constexpr std::uint32_t ControlFrames=32;
constexpr std::size_t TargetCount=8;
constexpr std::size_t WindowSize=1024;
// Catmull-Rom over a plain array, clamped at the ends.
[[nodiscard]] double readCubic(const std::vector<float>& t,double position) noexcept {
    const auto last=std::int64_t(t.size())-1;
    const double whole=std::floor(position);const double f=position-whole;
    const auto i=static_cast<std::int64_t>(whole);
    const auto at=[&](std::int64_t k){return double(t[std::size_t(std::clamp<std::int64_t>(k,0,last))]);};
    const double a=at(i-1),b=at(i),c=at(i+1),d=at(i+2);
    return b+.5*f*((c-a)+f*((2*a-5*b+4*c-d)+f*(3*(b-c)+d-a)));
}
struct Grain {double position{},step{},age{},length{},gain{};bool active{};};
struct MappingState {ParameterMapping mapping{};double value{},pole{};bool primed{};};
struct Layer {
    SampleLayer s;
    double rate{};                  // sample frames per output frame at pitch 1
    std::vector<MappingState> mappings;
    std::array<double,TargetCount> target{};
    // Loop region in sample frames
    double loopStart{},loopEnd{},fade{};
    double position{};
    bool playing{};
    double low{},lowPole{};
    double gainNow{},gainStep{};
    std::vector<Grain> grains;
    double untilGrain{};
    std::uint32_t liveGrains{};
};
}
double evaluateCurve(const ParameterCurve& c,double input) noexcept {
    const double span=c.inputMaximum-c.inputMinimum;
    double t=span!=0?(input-c.inputMinimum)/span:0;
    t=std::isfinite(t)?std::clamp(t,0.,1.):0;
    if(c.exponent!=1)t=std::pow(t,c.exponent);
    if(c.logarithmic&&c.outputMinimum>0&&c.outputMaximum>0)return c.outputMinimum*std::pow(c.outputMaximum/c.outputMinimum,t);
    return c.outputMinimum+(c.outputMaximum-c.outputMinimum)*t;
}
struct LayerInstrument::Impl {
    LayerInstrumentSettings s;
    std::vector<Layer> layers;
    std::array<std::atomic<double>,MaximumParameters> parameters{};
    std::atomic<std::uint32_t> triggers{0};
    std::atomic<double> releaseSeconds{-1};
    std::atomic<bool> done{false};
    std::uint32_t triggersSeen{};
    std::array<float,WindowSize+1> window{};
    std::uint64_t frame{},rng{};
    double fade{1},fadeStep{};
    bool started{};
    LayerStats stats;
    double uniform() noexcept {rng+=0x9e3779b97f4a7c15ull;std::uint64_t z=rng;z=(z^(z>>30))*0xbf58476d1ce4e5b9ull;z=(z^(z>>27))*0x94d049bb133111ebull;z^=z>>31;return (double(z>>11)+.5)/9007199254740992.;}
    void control(Layer& l) noexcept {
        auto& t=l.target;
        t={l.s.gain,l.s.pitch,l.s.lowPassHz,l.s.grainDensityPerSecond,l.s.grainSeconds,l.s.grainPosition,l.s.grainPositionJitter,l.s.grainPitchJitter};
        for(auto& m:l.mappings){
            const double wanted=evaluateCurve(m.mapping.curve,parameters[m.mapping.parameter].load(std::memory_order_relaxed));
            if(!m.primed){m.value=wanted;m.primed=true;}else m.value+=m.pole*(wanted-m.value);
            const auto index=std::size_t(m.mapping.target);
            if(m.mapping.target==LayerTarget::Gain||m.mapping.target==LayerTarget::Pitch)t[index]*=m.value;else t[index]=m.value;
        }
        const double cutoff=t[std::size_t(LayerTarget::LowPassHz)];
        l.lowPole=cutoff>0?1-std::exp(-2*Pi*std::min(cutoff,.45*s.sampleRateHz)/s.sampleRateHz):1;
        l.gainStep=(t[0]-l.gainNow)/ControlFrames;
    }
    void start() noexcept {
        started=true;fade=1;fadeStep=0;
        done.store(false,std::memory_order_relaxed);
        for(auto& l:layers){
            l.position=0;l.playing=true;l.untilGrain=0;
            for(auto& g:l.grains)g.active=false;
            l.liveGrains=0;
        }
    }
    void spawn(Layer& l) noexcept {
        const auto& t=l.target;
        const double density=t[std::size_t(LayerTarget::GrainDensity)];
        const double seconds=std::max(t[std::size_t(LayerTarget::GrainSeconds)],.001);
        // Poisson arrivals. Even spacing would be heard as a pitch.
        if(!(density>0)){l.untilGrain+=ControlFrames;return;}
        l.untilGrain+=-std::log(uniform())*s.sampleRateHz/density;
        Grain* g=nullptr;
        for(auto& candidate:l.grains)if(!candidate.active){g=&candidate;break;}
        if(!g){++stats.grainsRefused;return;}
        const double step=l.rate*t[1]*(1+t[std::size_t(LayerTarget::GrainPitchJitter)]*(2*uniform()-1));
        const double length=seconds*s.sampleRateHz;
        const double frames=double(l.s.sample->frames.size());
        const double room=std::max(frames-1-length*step,0.);
        const double centre=t[std::size_t(LayerTarget::GrainPosition)]+t[std::size_t(LayerTarget::GrainPositionJitter)]*(2*uniform()-1);
        // Uncorrelated grains add in power, and a Hann window keeps 3/8 of it.
        const double overlap=density*seconds*.375;
        *g={std::clamp(centre,0.,1.)*room,step,0,length,1/std::sqrt(std::max(overlap,1.)),true};
        ++l.liveGrains;++stats.grainsStarted;
    }
    void renderLayer(Layer& l,float* out,std::uint32_t frames) noexcept {
        const auto& t=l.s.sample->frames;
        double low=l.low,gain=l.gainNow;
        const double pole=l.lowPole,gainStep=l.gainStep,step=l.rate*l.target[1];
        for(std::uint32_t n=0;n<frames;++n){
            double value=0;
            if(l.s.kind==LayerKind::Granular){
                if(l.playing&&(l.untilGrain-=1)<=0)spawn(l);
                if(l.liveGrains)for(auto& g:l.grains){
                    if(!g.active)continue;
                    const double phase=g.age/g.length*WindowSize;
                    const auto index=static_cast<std::size_t>(phase);
                    const double w=double(window[index])+(double(window[index+1])-double(window[index]))*(phase-double(index));
                    value+=g.gain*w*readCubic(t,g.position);
                    g.position+=g.step;g.age+=1;
                    if(g.age>=g.length){g.active=false;--l.liveGrains;}
                }
            }else if(l.playing){
                value=readCubic(t,l.position);
                if(l.s.kind==LayerKind::Loop){
                    const double fadeStart=l.loopEnd-l.fade;
                    if(l.position>=fadeStart&&l.fade>0){
                        const double x=(l.position-fadeStart)/l.fade;
                        value=value*std::cos(.5*Pi*x)+readCubic(t,l.loopStart+(l.position-fadeStart))*std::sin(.5*Pi*x);
                    }
                    l.position+=step;
                    if(l.position>=l.loopEnd)l.position+=l.loopStart+l.fade-l.loopEnd;
                }else{
                    l.position+=step;
                    if(l.position>=double(t.size()))l.playing=false;
                }
            }
            low+=pole*(value-low);
            out[n]+=float(low*gain);
            gain+=gainStep;
        }
        l.low=low;l.gainNow=gain;
    }
};
LayerInstrument::LayerInstrument(std::vector<SampleLayer> layers,const LayerInstrumentSettings& s):impl_(std::make_unique<Impl>()){
    const auto bad=[](const std::string& what){throw std::invalid_argument("layer instrument: "+what);};
    if(!std::isfinite(s.sampleRateHz)||s.sampleRateHz<8000||s.sampleRateHz>384000)bad("sample rate must be in 8..384 kHz");
    if(s.maximumGrains<1||s.maximumGrains>1024)bad("maximum grains must be in 1..1024");
    if(layers.empty())bad("at least one layer is required");
    auto& m=*impl_;m.s=s;m.rng=(std::uint64_t(s.seed)<<32)|0x2545f491u;
    for(std::size_t i=0;i<=WindowSize;++i){const double x=std::sin(Pi*double(i)/WindowSize);m.window[i]=float(x*x);}
    for(auto& given:layers){
        if(!given.sample||given.sample->frames.size()<4||!(given.sample->sampleRateHz>0))bad("each layer needs a sample of at least 4 frames");
        if(!std::isfinite(given.pitch)||!(given.pitch>0)||!std::isfinite(given.gain)||given.gain<0)bad("gain or pitch invalid");
        Layer l;l.s=std::move(given);
        const auto& sample=*l.s.sample;
        l.rate=sample.sampleRateHz/s.sampleRateHz;
        for(const auto& mapping:l.s.mappings){
            if(mapping.parameter>=MaximumParameters)bad("parameter id out of range");
            if(!std::isfinite(mapping.smoothingSeconds)||mapping.smoothingSeconds<0)bad("smoothing must be non-negative");
            const double pole=mapping.smoothingSeconds>0?1-std::exp(-double(ControlFrames)/(mapping.smoothingSeconds*s.sampleRateHz)):1;
            l.mappings.push_back({mapping,0,pole,false});
        }
        if(l.s.kind==LayerKind::Loop){
            const double size=double(sample.frames.size());
            l.loopStart=sample.loop?double(sample.loop->startFrame):0;
            l.loopEnd=sample.loop?std::min(double(sample.loop->endFrame),size):size;
            if(!(l.loopEnd>l.loopStart+3))bad("loop region is too short");
            if(!std::isfinite(l.s.crossfadeSeconds)||l.s.crossfadeSeconds<0)bad("crossfade must be non-negative");
            l.fade=std::min(l.s.crossfadeSeconds*sample.sampleRateHz,(l.loopEnd-l.loopStart)*.5);
        }
        if(l.s.kind==LayerKind::Granular)l.grains.assign(s.maximumGrains,Grain{});
        m.layers.push_back(std::move(l));
    }
}
LayerInstrument::~LayerInstrument()=default;
LayerInstrument::LayerInstrument(LayerInstrument&&) noexcept=default;
LayerInstrument& LayerInstrument::operator=(LayerInstrument&&) noexcept=default;
void LayerInstrument::setParameter(ParameterId p,double value) noexcept {if(p<MaximumParameters&&std::isfinite(value))impl_->parameters[p].store(value,std::memory_order_relaxed);}
void LayerInstrument::trigger() noexcept {impl_->done.store(false,std::memory_order_relaxed);impl_->triggers.fetch_add(1,std::memory_order_release);}
void LayerInstrument::release(double seconds) noexcept {impl_->releaseSeconds.store(seconds>0?seconds:0,std::memory_order_release);}
void LayerInstrument::render(std::span<float> mono) noexcept {
    auto& m=*impl_;
    std::fill(mono.begin(),mono.end(),0.f);
    std::size_t done=0;
    while(done<mono.size()){
        const auto into=std::uint32_t(m.frame%ControlFrames);
        const auto chunk=static_cast<std::uint32_t>(std::min<std::size_t>(ControlFrames-into,mono.size()-done));
        if(into==0){
            if(const auto t=m.triggers.load(std::memory_order_acquire);t!=m.triggersSeen){m.triggersSeen=t;m.start();}
            if(const auto r=m.releaseSeconds.exchange(-1,std::memory_order_acq_rel);r>=0)m.fadeStep=r>0?-1/(r*m.s.sampleRateHz):-1;
            for(auto& l:m.layers)m.control(l);
        }
        for(auto& l:m.layers)m.renderLayer(l,mono.data()+done,chunk);
        for(std::uint32_t n=0;n<chunk;++n){
            mono[done+n]*=float(m.fade);
            if(m.fadeStep!=0)m.fade=std::max(m.fade+m.fadeStep,0.);
        }
        if(m.fade==0)for(auto& l:m.layers){
            l.playing=false;l.liveGrains=0;
            for(auto& g:l.grains)g.active=false;
        }
        done+=chunk;m.frame+=chunk;
    }
    if(m.started){
        bool live=false;
        for(const auto& l:m.layers)live=live||l.playing||l.liveGrains;
        if(!live&&m.triggers.load(std::memory_order_relaxed)==m.triggersSeen)m.done.store(true,std::memory_order_relaxed);
    }
}
bool LayerInstrument::finished() const noexcept {return impl_->done.load(std::memory_order_relaxed);}
const LayerStats& LayerInstrument::stats() const noexcept {return impl_->stats;}
}
