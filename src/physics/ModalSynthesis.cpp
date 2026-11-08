// Runtime modal synthesis: per-mode two-pole recursive damped-sinusoid
// resonators in double state, excited by half-sine contact pulses or a
// continuous contact force, culled per voice, with energy-based voice stealing.
#include "marathon_spatial/physics/Modal.hpp"
#include <algorithm>
#include <atomic>
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <string>

namespace mspatial::physics {
namespace {[[nodiscard]] bool finite(Vec3 v){return std::isfinite(v.x)&&std::isfinite(v.y)&&std::isfinite(v.z);}}
double halfSinePulseSpectrum(double omega,double tau){
    if(!std::isfinite(omega)||omega<0)throw std::invalid_argument("half-sine spectrum: angular frequency must be finite and non-negative");
    if(!std::isfinite(tau)||!(tau>0))throw std::invalid_argument("half-sine spectrum: contact duration must be finite and positive");
    const double x=omega*tau/std::numbers::pi,e=x-1;
    // Removable singularity at x = 1: cos(pi x/2)/(1-x^2) = sin(pi e/2)/(e(2+e)).
    if(std::fabs(e)<1e-3){const double h=std::numbers::pi*e/2;return (std::numbers::pi/2)*(1-h*h/6)/(2+e);}
    return std::cos(std::numbers::pi*x/2)/(1-x*x);
}
double hertzContactSeconds(const HertzContact& c){
    const auto bad=[](const char* what){throw std::invalid_argument(std::string("Hertz contact: ")+what);};
    if(!std::isfinite(c.relativeNormalSpeed)||!(c.relativeNormalSpeed>0))bad("relative normal speed must be finite and positive");
    if(!std::isfinite(c.massAKg)||!(c.massAKg>0))bad("mass A must be finite and positive");
    if(std::isnan(c.massBKg)||!(c.massBKg>0))bad("mass B must be positive (infinity for a static body)");
    if(!std::isfinite(c.effectiveRadiusMeters)||!(c.effectiveRadiusMeters>0))bad("effective radius must be finite and positive");
    for(double e:{c.youngModulusAPa,c.youngModulusBPa})if(!std::isfinite(e)||!(e>0))bad("Young modulus must be finite and positive");
    for(double nu:{c.poissonRatioA,c.poissonRatioB})if(!std::isfinite(nu)||nu<0||!(nu<.5))bad("Poisson ratio must be in [0, 0.5)");
    const double mass=std::isinf(c.massBKg)?c.massAKg:c.massAKg*c.massBKg/(c.massAKg+c.massBKg);
    const double modulus=1/((1-c.poissonRatioA*c.poissonRatioA)/c.youngModulusAPa+(1-c.poissonRatioB*c.poissonRatioB)/c.youngModulusBPa);
    const double tau=2.87*std::pow(mass*mass/(c.effectiveRadiusMeters*modulus*modulus*c.relativeNormalSpeed),.2);
    return std::clamp(tau,MinimumContactSeconds,MaximumContactSeconds);
}

namespace {
constexpr std::uint16_t NoSlot=0xFFFF;
constexpr std::uint32_t Lanes=8; // render kernel width; inactive padding slots hold zeros
struct PreparedMode{
    double omega{},omegaDamped{},r{},sinTheta{},cosTheta{},c1{},c2{},injectPrevious{},injectBeforePrevious{},weight{},forceToState{};
    bool renderable{};
};
struct PreparedModel{std::shared_ptr<const ModalModel> model;std::vector<PreparedMode> modes;};
// One voice's continuous force. Gains ramp per frame, so a Continue cadence of
// 10 Hz gives no zipper; `rng` is seeded from the instance (deterministic).
struct DriveState{
    double gain{},gainTarget{},gainStep{},pole{},poleTarget{},poleStep{},low{};
    std::uint64_t rng{};
    const float* external{};  // this render's caller-owned force, N per frame
    bool active{},usesExternal{};
};
[[nodiscard]] double unitNoise(std::uint64_t& state){
    state+=0x9e3779b97f4a7c15ull;std::uint64_t z=state;z=(z^(z>>30))*0xbf58476d1ce4e5b9ull;z=(z^(z>>27))*0x94d049bb133111ebull;z^=z>>31;
    return (double(z>>11)*(1./9007199254740992.)*2-1)*1.7320508075688772; // uniform, unit variance
}
}

struct ModalSynth::Impl{
    ModalSynthSettings s;
    std::uint32_t capacity{}; // per-voice slots padded to a multiple of Lanes
    // Fixed at construction: registration writes the next entry, then publishes the count, so a reader on another thread sees only
    // complete models and none ever moves.
    std::unique_ptr<PreparedModel[]> models;
    std::atomic<std::uint32_t> modelCount{0};
    [[nodiscard]] std::uint32_t modelsPublished() const noexcept {return modelCount.load(std::memory_order_acquire);}
    std::vector<std::uint8_t> active;
    std::vector<std::uint64_t> instance;
    std::vector<ModalModelHandle> model;
    std::vector<std::uint32_t> count;
    std::vector<double> c1,c2,y1,y2,w,in;  // voices x capacity
    std::vector<double> perceptual;        // per voice
    std::vector<DriveState> drives;        // per voice
    std::vector<std::uint16_t> modeOf;     // voices x capacity
    std::vector<std::uint16_t> slotOf;     // voices x MaximumModes
    std::vector<double> candidateAmplitude,candidateOutput;std::vector<std::uint16_t> order;
    std::vector<double> accumulator,busAccumulator,force;

    [[nodiscard]] std::size_t at(std::uint32_t v,std::uint32_t j) const {return std::size_t(v)*capacity+j;}
    [[nodiscard]] double slotAmplitude(std::uint32_t v,std::uint32_t j) const {
        const auto& p=models[model[v]].modes[modeOf[at(v,j)]];const double a=y1[at(v,j)],b=y2[at(v,j)];
        const double q=(a*p.cosTheta-p.r*b)/p.sinTheta;
        return std::fabs(w[at(v,j)])*std::sqrt(a*a+q*q);
    }
    [[nodiscard]] double energy(std::uint32_t v) const {
        double e=0;
        for(std::uint32_t j=0;j<count[v];++j){const double a=slotAmplitude(v,j);e+=a*a/2;}
        return e;
    }
    // What a voice defends when voices compete. Driven voices hold their slot.
    [[nodiscard]] double defended(std::uint32_t v) const {
        return drives[v].active&&drives[v].gainTarget>0?std::numeric_limits<double>::infinity():energy(v)*perceptual[v]*perceptual[v];
    }
    void clearSlot(std::uint32_t v,std::uint32_t j){const auto i=at(v,j);c1[i]=c2[i]=y1[i]=y2[i]=w[i]=in[i]=0;modeOf[i]=0;}
    void reset(std::uint32_t v){
        for(std::uint32_t j=0;j<count[v];++j){slotOf[std::size_t(v)*ModalAnalysisSettings::MaximumModes+modeOf[at(v,j)]]=NoSlot;clearSlot(v,j);}
        count[v]=0;active[v]=0;instance[v]=0;model[v]=0;perceptual[v]=1;drives[v]={};
    }
    // The voice of an instance on a model: its own, a free one, or the quietest
    // when `incoming` beats it. Returns maximumVoices when refused.
    [[nodiscard]] std::uint32_t acquire(std::uint64_t key,ModalModelHandle handle,double incoming,double weight,ModalExciteResult& result){
        std::uint32_t v=s.maximumVoices;
        for(std::uint32_t i=0;i<s.maximumVoices;++i)if(active[i]&&instance[i]==key){v=i;break;}
        if(v<s.maximumVoices&&model[v]==handle){result.outcome=ModalExciteOutcome::Merged;perceptual[v]=std::max(perceptual[v],weight);return v;}
        if(v<s.maximumVoices)reset(v); // the instance now sounds a different model
        else{
            for(std::uint32_t i=0;i<s.maximumVoices;++i)if(!active[i]){v=i;break;}
            if(v==s.maximumVoices){
                double quietest=std::numeric_limits<double>::infinity();
                for(std::uint32_t i=0;i<s.maximumVoices;++i){const double en=defended(i);if(en<quietest){quietest=en;v=i;}}
                if(v==s.maximumVoices||!(incoming*weight*weight>quietest)){result.outcome=ModalExciteOutcome::RejectedQuieter;result.excitationEnergy=incoming;return s.maximumVoices;}
                result.outcome=ModalExciteOutcome::Stole;result.stolenInstance=instance[v];reset(v);
            }
        }
        if(result.outcome!=ModalExciteOutcome::Stole)result.outcome=ModalExciteOutcome::NewVoice;
        active[v]=1;instance[v]=key;model[v]=handle;count[v]=0;perceptual[v]=weight;
        return v;
    }
    void removeSlot(std::uint32_t v,std::uint32_t j){
        const std::uint32_t last=count[v]-1;
        slotOf[std::size_t(v)*ModalAnalysisSettings::MaximumModes+modeOf[at(v,j)]]=NoSlot;
        if(j!=last){
            const auto from=at(v,last),to=at(v,j);
            c1[to]=c1[from];c2[to]=c2[from];y1[to]=y1[from];y2[to]=y2[from];w[to]=w[from];in[to]=in[from];modeOf[to]=modeOf[from];
            slotOf[std::size_t(v)*ModalAnalysisSettings::MaximumModes+modeOf[to]]=static_cast<std::uint16_t>(j);
        }
        clearSlot(v,last);--count[v];
    }
};

ModalSynth::ModalSynth(const ModalSynthSettings& s):impl_(std::make_unique<Impl>()){
    const auto bad=[](const std::string& what){throw std::invalid_argument("modal synth: "+what);};
    if(!std::isfinite(s.sampleRateHz)||s.sampleRateHz<8000||s.sampleRateHz>384000)bad("sample rate must be in 8..384 kHz");
    if(s.maximumVoices<1||s.maximumVoices>ModalSynthSettings::MaximumVoices)bad("maximum voices must be in 1.."+std::to_string(ModalSynthSettings::MaximumVoices));
    if(s.maximumModesPerVoice<1||s.maximumModesPerVoice>ModalSynthSettings::MaximumModesPerVoice)bad("maximum modes per voice must be in 1.."+std::to_string(ModalSynthSettings::MaximumModesPerVoice));
    if(s.maximumModels<1||s.maximumModels>ModalSynthSettings::MaximumModels)bad("maximum models must be in 1.."+std::to_string(ModalSynthSettings::MaximumModels));
    if(s.maximumBlockFrames<1||s.maximumBlockFrames>ModalSynthSettings::MaximumBlockFrames)bad("maximum block frames must be in 1.."+std::to_string(ModalSynthSettings::MaximumBlockFrames));
    if(!std::isfinite(s.cullAmplitude)||s.cullAmplitude<0)bad("cull amplitude must be finite and non-negative");
    if(!std::isfinite(s.outputGain)||!(s.outputGain>0)||s.outputGain>1e6)bad("output gain must be finite in (0, 1e6]");
    auto& m=*impl_;m.s=s;m.capacity=(s.maximumModesPerVoice+Lanes-1)/Lanes*Lanes;
    m.models=std::make_unique<PreparedModel[]>(s.maximumModels);
    const std::size_t slots=std::size_t(s.maximumVoices)*m.capacity;
    m.active.assign(s.maximumVoices,0);m.instance.assign(s.maximumVoices,0);m.model.assign(s.maximumVoices,0);m.count.assign(s.maximumVoices,0);
    for(auto* v:{&m.c1,&m.c2,&m.y1,&m.y2,&m.w,&m.in})v->assign(slots,0.);
    m.perceptual.assign(s.maximumVoices,1.);m.drives.assign(s.maximumVoices,DriveState{});
    m.modeOf.assign(slots,0);m.slotOf.assign(std::size_t(s.maximumVoices)*ModalAnalysisSettings::MaximumModes,NoSlot);
    m.candidateAmplitude.assign(ModalAnalysisSettings::MaximumModes,0.);m.candidateOutput.assign(ModalAnalysisSettings::MaximumModes,0.);m.order.assign(ModalAnalysisSettings::MaximumModes,0);
    m.accumulator.assign(s.maximumBlockFrames,0.);m.busAccumulator.assign(s.maximumBlockFrames,0.);m.force.assign(s.maximumBlockFrames,0.);
}
ModalSynth::~ModalSynth()=default;
ModalSynth::ModalSynth(ModalSynth&&) noexcept=default;
ModalSynth& ModalSynth::operator=(ModalSynth&&) noexcept=default;
const ModalSynthSettings& ModalSynth::settings() const noexcept {return impl_->s;}

ModalModelHandle ModalSynth::registerModel(std::shared_ptr<const ModalModel> model){
    if(!model)throw std::invalid_argument("modal synth: null model");
    validateModalModel(*model);
    auto& m=*impl_;
    const std::uint32_t next=m.modelCount.load(std::memory_order_relaxed);
    if(next>=m.s.maximumModels)throw std::invalid_argument("modal synth: model table full");
    PreparedModel p;p.modes.resize(model->modes.size());
    const double T=1/m.s.sampleRateHz;
    for(std::size_t k=0;k<model->modes.size();++k){
        const auto& mode=model->modes[k];auto& q=p.modes[k];
        q.omega=mode.angularFrequency;q.omegaDamped=mode.angularFrequency*std::sqrt(1-mode.dampingRatio*mode.dampingRatio);
        const double theta=q.omegaDamped*T;
        // Modes above 0.45 fs cannot be rendered without aliasing at this rate.
        q.renderable=q.omegaDamped>0&&q.omegaDamped/(2*std::numbers::pi)<=.45*m.s.sampleRateHz;
        if(!q.renderable)continue;
        q.r=std::exp(-mode.dampingRatio*mode.angularFrequency*T);q.sinTheta=std::sin(theta);q.cosTheta=std::cos(theta);
        q.c1=2*q.r*q.cosTheta;q.c2=q.r*q.r;
        // State (q[-1], q[-2]) of A r^n sin(theta n) per unit A: the sequence
        // starts at q[0] = 0 on the next rendered frame.
        q.injectPrevious=-q.sinTheta/q.r;q.injectBeforePrevious=-std::sin(2*theta)/(q.r*q.r);
        // Output: surface-RMS normal velocity proxy, outputGain x R_k x omega_k x q_k.
        q.weight=m.s.outputGain*mode.surfaceRms*mode.angularFrequency;
        // One frame of force F is the impulse F T: state input T r sin(theta) / omega_d per unit gain.
        q.forceToState=T*q.r*q.sinTheta/q.omegaDamped;
    }
    p.model=std::move(model);
    m.models[next]=std::move(p);
    m.modelCount.store(next+1,std::memory_order_release);
    return static_cast<ModalModelHandle>(next);
}

ModalExciteResult ModalSynth::excite(const ModalExcitation& e){
    auto& m=*impl_;
    const auto bad=[](const std::string& what){throw std::invalid_argument("modal excite: "+what);};
    if(e.instance==0)bad("instance key must be nonzero");
    if(e.model>=m.modelsPublished())bad("unknown model handle");
    if(!std::isfinite(e.impulseNs)||e.impulseNs<0)bad("impulse must be finite and non-negative");
    if(!std::isfinite(e.contactSeconds)||e.contactSeconds<MinimumContactSeconds||e.contactSeconds>MaximumContactSeconds)bad("contact duration outside 1 us..100 ms");
    if(!std::isfinite(e.perceptualWeight)||!(e.perceptualWeight>0))bad("perceptual weight must be finite and positive");
    const auto& prepared=m.models[e.model];const ModalModel& model=*prepared.model;
    ModalExciteResult result;
    if(e.sample){if(*e.sample>=model.samples.size())bad("contact sample out of range");result.sample=*e.sample;}
    else{if(!finite(e.localPosition))bad("contact position is not finite");result.sample=nearestModalSample(model,e.localPosition);}
    // Per-mode amplitude = gain x impulse x pulse spectrum / (omega_d m), modal mass 1.
    std::uint32_t candidates=0;
    const double floor=m.s.cullAmplitude;
    for(std::size_t k=0;k<model.modes.size();++k){
        const auto& p=prepared.modes[k];
        if(!p.renderable)continue;
        const double a=double(model.gain(result.sample,k))*e.impulseNs*halfSinePulseSpectrum(p.omega,e.contactSeconds)/p.omegaDamped;
        const double out=std::fabs(p.weight*a);
        if(!(out>floor))continue;
        m.candidateAmplitude[k]=a;m.candidateOutput[k]=out;m.order[candidates++]=static_cast<std::uint16_t>(k);
    }
    if(candidates==0)return result;
    std::sort(m.order.begin(),m.order.begin()+candidates,[&](std::uint16_t a,std::uint16_t b){
        return m.candidateOutput[a]!=m.candidateOutput[b]?m.candidateOutput[a]>m.candidateOutput[b]:a<b;});
    if(e.maximumModes)candidates=std::min(candidates,e.maximumModes);
    const std::uint32_t limit=std::min(candidates,m.s.maximumModesPerVoice);
    double incoming=0;for(std::uint32_t i=0;i<limit;++i)incoming+=m.candidateOutput[m.order[i]]*m.candidateOutput[m.order[i]]/2;

    const std::uint32_t v=m.acquire(e.instance,e.model,incoming,e.perceptualWeight,result);
    if(v==m.s.maximumVoices)return result;
    result.voice=v;
    std::uint16_t* slotOf=m.slotOf.data()+std::size_t(v)*ModalAnalysisSettings::MaximumModes;
    bool saturated=false; // a full voice refused a candidate; weaker new modes are refused too
    for(std::uint32_t i=0;i<candidates;++i){
        const std::uint16_t k=m.order[i];const auto& p=prepared.modes[k];
        std::uint32_t j=slotOf[k];
        if(j==NoSlot){
            if(saturated)continue; // modes already sounding still receive their share (linearity)
            if(m.count[v]<m.s.maximumModesPerVoice)j=m.count[v]++;
            else{
                std::uint32_t weakest=0;double weakestAmplitude=std::numeric_limits<double>::infinity();
                for(std::uint32_t q=0;q<m.count[v];++q){const double a=m.slotAmplitude(v,q);if(a<weakestAmplitude){weakestAmplitude=a;weakest=q;}}
                if(!(m.candidateOutput[k]>weakestAmplitude)){saturated=true;continue;}
                slotOf[m.modeOf[m.at(v,weakest)]]=NoSlot;m.clearSlot(v,weakest);j=weakest;
            }
            const auto s=m.at(v,j);
            m.c1[s]=p.c1;m.c2[s]=p.c2;m.w[s]=p.weight;m.modeOf[s]=k;slotOf[k]=static_cast<std::uint16_t>(j);
        }
        const auto s=m.at(v,j);const double a=m.candidateAmplitude[k];
        m.y1[s]+=a*p.injectPrevious;m.y2[s]+=a*p.injectBeforePrevious;
        ++result.modesExcited;result.excitationEnergy+=m.candidateOutput[k]*m.candidateOutput[k]/2;
    }
    return result;
}

ModalExciteResult ModalSynth::drive(const ModalDrive& d){
    auto& m=*impl_;
    const auto bad=[](const std::string& what){throw std::invalid_argument("modal drive: "+what);};
    if(d.instance==0)bad("instance key must be nonzero");
    if(d.model>=m.modelsPublished())bad("unknown model handle");
    if(!std::isfinite(d.forceRmsN)||d.forceRmsN<0)bad("force must be finite and non-negative");
    if(!std::isfinite(d.cutoffHz)||!(d.cutoffHz>0))bad("cutoff must be finite and positive");
    if(!std::isfinite(d.rampSeconds)||d.rampSeconds<0)bad("ramp must be finite and non-negative");
    if(!std::isfinite(d.perceptualWeight)||!(d.perceptualWeight>0))bad("perceptual weight must be finite and positive");
    const auto& prepared=m.models[d.model];const ModalModel& model=*prepared.model;
    ModalExciteResult result;
    if(d.sample){if(*d.sample>=model.samples.size())bad("contact sample out of range");result.sample=*d.sample;}
    else{if(!finite(d.localPosition))bad("contact position is not finite");result.sample=nearestModalSample(model,d.localPosition);}
    const double fs=m.s.sampleRateHz;
    if(d.forceRmsN==0){ // Release: a voice that was never driven is left alone.
        for(std::uint32_t v=0;v<m.s.maximumVoices;++v)if(m.active[v]&&m.instance[v]==d.instance&&m.drives[v].active){
            auto& state=m.drives[v];const double frames=std::max(1.,d.rampSeconds*fs);
            state.gainTarget=0;state.gainStep=(0-state.gain)/frames;result.outcome=ModalExciteOutcome::Merged;result.voice=v;}
        return result;
    }
    // Rank by the output one unit of force reaches: the drive's share of a mode is its gain at the sample.
    std::uint32_t candidates=0;
    for(std::size_t k=0;k<model.modes.size();++k){
        const auto& p=prepared.modes[k];
        if(!p.renderable)continue;
        const double response=std::fabs(p.weight*double(model.gain(result.sample,k))*p.forceToState);
        if(!(response>0))continue;
        m.candidateOutput[k]=response;m.order[candidates++]=static_cast<std::uint16_t>(k);
    }
    if(candidates==0)return result;
    std::sort(m.order.begin(),m.order.begin()+candidates,[&](std::uint16_t a,std::uint16_t b){
        return m.candidateOutput[a]!=m.candidateOutput[b]?m.candidateOutput[a]>m.candidateOutput[b]:a<b;});
    candidates=std::min({candidates,m.s.maximumModesPerVoice,d.maximumModes?d.maximumModes:m.s.maximumModesPerVoice});
    double incoming=0;for(std::uint32_t i=0;i<candidates;++i){const double a=m.candidateOutput[m.order[i]]*d.forceRmsN;incoming+=a*a/2;}
    const std::uint32_t v=m.acquire(d.instance,d.model,incoming,d.perceptualWeight,result);
    if(v==m.s.maximumVoices)return result;
    result.voice=v;
    std::uint16_t* slotOf=m.slotOf.data()+std::size_t(v)*ModalAnalysisSettings::MaximumModes;
    for(std::uint32_t j=0;j<m.count[v];++j)m.in[m.at(v,j)]=0; // the contact may have moved to another sample
    for(std::uint32_t i=0;i<candidates;++i){
        const std::uint16_t k=m.order[i];const auto& p=prepared.modes[k];
        std::uint32_t j=slotOf[k];
        if(j==NoSlot){
            if(m.count[v]<m.s.maximumModesPerVoice)j=m.count[v]++;
            else{ // Displace the weakest ringing mode this contact does not feed.
                std::uint32_t weakest=m.count[v];double weakestAmplitude=std::numeric_limits<double>::infinity();
                for(std::uint32_t q=0;q<m.count[v];++q){if(m.in[m.at(v,q)]!=0)continue;const double a=m.slotAmplitude(v,q);if(a<weakestAmplitude){weakestAmplitude=a;weakest=q;}}
                if(weakest==m.count[v])continue;
                slotOf[m.modeOf[m.at(v,weakest)]]=NoSlot;m.clearSlot(v,weakest);j=weakest;
            }
            const auto s=m.at(v,j);
            m.c1[s]=p.c1;m.c2[s]=p.c2;m.w[s]=p.weight;m.modeOf[s]=k;slotOf[k]=static_cast<std::uint16_t>(j);
        }
        m.in[m.at(v,j)]=double(model.gain(result.sample,k))*p.forceToState;
        ++result.modesExcited;
    }
    result.excitationEnergy=incoming;
    auto& state=m.drives[v];
    const double frames=std::max(1.,d.rampSeconds*fs);
    const double pole=1-std::exp(-2*std::numbers::pi*std::min(d.cutoffHz,.45*fs)/fs);
    if(!state.active){state.pole=pole;state.low=0;state.rng=d.instance*0x9e3779b97f4a7c15ull+1;state.active=true;}
    state.usesExternal=d.externalForce;
    const double target=d.externalForce?1.:d.forceRmsN; // an external force carries its own level; the gain is its envelope
    state.gainTarget=target;state.gainStep=(target-state.gain)/frames;
    state.poleTarget=pole;state.poleStep=(pole-state.pole)/frames;
    return result;
}

void ModalSynth::setDriveForce(std::uint32_t slot,const float* force) noexcept {
    auto& m=*impl_;
    if(slot<m.s.maximumVoices&&m.drives[slot].active)m.drives[slot].external=force;
}
bool ModalSynth::driven(std::uint32_t slot) const noexcept {
    const auto& m=*impl_;
    return slot<m.s.maximumVoices&&m.active[slot]&&m.drives[slot].active;
}
std::uint32_t ModalSynth::modelCount() const noexcept {return impl_->modelsPublished();}
std::uint32_t ModalSynth::sampleCount(ModalModelHandle model) const noexcept {
    const auto& m=*impl_;
    return model<m.modelsPublished()?static_cast<std::uint32_t>(m.models[model].model->samples.size()):0;
}
std::uint64_t ModalSynth::instanceIn(std::uint32_t slot) const noexcept {
    const auto& m=*impl_;
    return slot<m.s.maximumVoices&&m.active[slot]?m.instance[slot]:0;
}
void ModalSynth::stop(std::uint64_t instance) noexcept {
    auto& m=*impl_;
    for(std::uint32_t v=0;v<m.s.maximumVoices;++v)if(m.active[v]&&m.instance[v]==instance)m.reset(v);
}

void ModalSynth::render(std::span<float> bus,std::span<float> voices,std::size_t frames){
    auto& m=*impl_;
    if(frames>m.s.maximumBlockFrames)throw std::invalid_argument("modal render: frames exceed maximumBlockFrames");
    if(bus.size()!=frames)throw std::invalid_argument("modal render: bus size must equal frames");
    if(!voices.empty()&&voices.size()!=std::size_t(m.s.maximumVoices)*frames)throw std::invalid_argument("modal render: voice buffer must be empty or maximumVoices x frames");
    double* acc=m.accumulator.data();double* busAcc=m.busAccumulator.data();
    std::fill_n(busAcc,frames,0.);
    std::fill(voices.begin(),voices.end(),0.f);
    for(std::uint32_t v=0;v<m.s.maximumVoices;++v){
        if(!m.active[v])continue;
        std::fill_n(acc,frames,0.);
        // Eight independent resonators per pass hide the multiply-subtract
        // latency; the per-frame sum has a fixed tree order (deterministic).
        const std::uint32_t blocks=(m.count[v]+Lanes-1)/Lanes;
        auto& drive=m.drives[v];
        if(drive.active){
            // Normalized one-pole noise: the force keeps its RMS at any cutoff.
            double* f=m.force.data();
            const float* external=drive.usesExternal?drive.external:nullptr;
            drive.external=nullptr;
            for(std::size_t n=0;n<frames;++n){
                if(drive.gain!=drive.gainTarget){drive.gain+=drive.gainStep;if((drive.gainStep>0&&drive.gain>drive.gainTarget)||(drive.gainStep<0&&drive.gain<drive.gainTarget)||drive.gainStep==0)drive.gain=drive.gainTarget;}
                if(drive.pole!=drive.poleTarget){drive.pole+=drive.poleStep;if((drive.poleStep>0&&drive.pole>drive.poleTarget)||(drive.poleStep<0&&drive.pole<drive.poleTarget)||drive.poleStep==0)drive.pole=drive.poleTarget;}
                if(drive.usesExternal){f[n]=external?drive.gain*double(external[n]):0.;continue;}
                drive.low+=drive.pole*(unitNoise(drive.rng)-drive.low);
                f[n]=drive.gain*drive.low*std::sqrt((2-drive.pole)/drive.pole);
            }
            if(drive.gain==0&&drive.gainTarget==0){drive=DriveState{};for(std::uint32_t j=0;j<m.count[v];++j)m.in[m.at(v,j)]=0;}
        }
        const bool driven=drive.active;
        for(std::uint32_t b=0;b<blocks;++b){
            const std::size_t base=m.at(v,Lanes*b);
            double c1[Lanes],c2[Lanes],w[Lanes],p[Lanes],q[Lanes],in[Lanes];
            for(std::uint32_t j=0;j<Lanes;++j){c1[j]=m.c1[base+j];c2[j]=m.c2[base+j];w[j]=m.w[base+j];p[j]=m.y1[base+j];q[j]=m.y2[base+j];in[j]=m.in[base+j];}
            const double* f=m.force.data();
            for(std::size_t n=0;n<frames;++n){
                double y[Lanes];
                if(driven)for(std::uint32_t j=0;j<Lanes;++j)y[j]=c1[j]*p[j]-c2[j]*q[j]+in[j]*f[n];
                else for(std::uint32_t j=0;j<Lanes;++j)y[j]=c1[j]*p[j]-c2[j]*q[j];
                acc[n]+=((w[0]*y[0]+w[1]*y[1])+(w[2]*y[2]+w[3]*y[3]))+((w[4]*y[4]+w[5]*y[5])+(w[6]*y[6]+w[7]*y[7]));
                for(std::uint32_t j=0;j<Lanes;++j){q[j]=p[j];p[j]=y[j];}
            }
            for(std::uint32_t j=0;j<Lanes;++j){m.y1[base+j]=p[j];m.y2[base+j]=q[j];}
        }
        float* out=voices.empty()?nullptr:voices.data()+std::size_t(v)*frames;
        for(std::size_t n=0;n<frames;++n){busAcc[n]+=acc[n];if(out)out[n]=static_cast<float>(acc[n]);}
        if(m.s.cullAmplitude>0&&!driven){
            for(std::uint32_t j=0;j<m.count[v];)if(m.slotAmplitude(v,j)<m.s.cullAmplitude)m.removeSlot(v,j);else ++j;
            if(m.count[v]==0)m.reset(v);
        }
    }
    for(std::size_t n=0;n<frames;++n)bus[n]=static_cast<float>(busAcc[n]);
}

ModalVoiceState ModalSynth::voice(std::uint32_t slot) const {
    const auto& m=*impl_;
    if(slot>=m.s.maximumVoices)throw std::invalid_argument("modal synth: voice slot out of range");
    if(!m.active[slot])return {};
    return {true,m.instance[slot],m.model[slot],m.count[slot],m.energy(slot),m.drives[slot].active};
}
std::uint32_t ModalSynth::activeVoices() const noexcept {
    std::uint32_t n=0;for(auto a:impl_->active)n+=a?1u:0u;return n;
}
}
