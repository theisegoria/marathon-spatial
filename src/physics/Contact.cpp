// Contact forces and the sample-accurate contact voice manager.
#include "marathon_spatial/physics/Contact.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <string>

namespace mspatial::physics {
namespace {
constexpr double Pi=std::numbers::pi;
[[nodiscard]] double unitGaussian(std::uint64_t& state) noexcept {
    const auto next=[&]{state+=0x9e3779b97f4a7c15ull;std::uint64_t z=state;z=(z^(z>>30))*0xbf58476d1ce4e5b9ull;z=(z^(z>>27))*0x94d049bb133111ebull;z^=z>>31;return (double(z>>11)+.5)/9007199254740992.;};
    const double u=next(),v=next();
    return std::sqrt(-2*std::log(u))*std::cos(2*Pi*v);
}
[[nodiscard]] std::uint32_t seedOf(std::uint64_t key) noexcept {return static_cast<std::uint32_t>((key*0x9e3779b97f4a7c15ull)>>32);}
[[nodiscard]] float lerpWrapped(const std::vector<float>& t,double position) noexcept {
    const std::size_t n=t.size();
    const double i=std::floor(position);const double f=position-i;
    const auto a=static_cast<std::size_t>(static_cast<std::int64_t>(i))&(n-1);
    return float(double(t[a])+(double(t[(a+1)&(n-1)])-double(t[a]))*f);
}
}

SurfaceProfile::SurfaceProfile(double cornerWavelengthM,std::uint32_t seed,std::uint32_t sizeLog2):wavelength_(cornerWavelengthM){
    if(!std::isfinite(cornerWavelengthM)||!(cornerWavelengthM>0))throw std::invalid_argument("surface profile: corner wavelength must be finite and positive");
    if(sizeLog2<10||sizeLog2>22)throw std::invalid_argument("surface profile: size outside 2^10..2^22");
    spacing_=cornerWavelengthM/32;
    const std::size_t n=std::size_t(1)<<sizeLog2;
    // AR(1) along the surface, closed into a loop by blending the last stretch into the first.
    const double rho=std::exp(-2*Pi*spacing_/cornerWavelengthM),drive=std::sqrt(1-rho*rho);
    const std::size_t blend=std::min<std::size_t>(n/8,256);
    std::vector<double> h(n+blend);
    std::uint64_t rng=(std::uint64_t(seed)<<32)|0x5851f42du;
    double value=unitGaussian(rng);
    for(auto& v:h){v=value;value=rho*value+drive*unitGaussian(rng);}
    for(std::size_t i=0;i<blend;++i){
        const double t=(i+.5)/double(blend),w=std::sin(.5*Pi*t);
        h[i]=std::cos(.5*Pi*t)*h[n+i]+w*h[i];
    }
    h.resize(n);
    double mean=0,square=0;
    for(double v:h)mean+=v;
    mean/=double(n);
    for(double& v:h){v-=mean;square+=v*v;}
    const double scale=1/std::sqrt(square/double(n));
    std::vector<float> base(n);
    for(std::size_t i=0;i<n;++i)base[i]=float(h[i]*scale);
    levels_.push_back(std::move(base));
    while(levels_.back().size()>16){
        const auto& up=levels_.back();
        std::vector<float> down(up.size()/2);
        for(std::size_t i=0;i<down.size();++i)down[i]=(up[2*i]+up[2*i+1])*.5f;
        levels_.push_back(std::move(down));
    }
}
float SurfaceProfile::read(double x,double width) const noexcept {
    const double position=x/spacing_;
    const double lod=width>spacing_?std::log2(width/spacing_):0;
    const double top=double(levels_.size()-1);
    if(lod>=top)return lerpWrapped(levels_.back(),position/std::ldexp(1.,int(top)));
    const auto low=static_cast<int>(lod);const double f=lod-low;
    const double a=lerpWrapped(levels_[std::size_t(low)],std::ldexp(position,-low));
    if(f==0)return float(a);
    const double b=lerpWrapped(levels_[std::size_t(low)+1],std::ldexp(position,-low-1));
    return float(a+(b-a)*f);
}

double effectiveModulusPa(const ModalElastic& a,const ModalElastic& b){
    if(!(a.youngModulusPa>0)||!(b.youngModulusPa>0))throw std::invalid_argument("effective modulus: Young modulus must be positive");
    return 1/((1-a.poissonRatio*a.poissonRatio)/a.youngModulusPa+(1-b.poissonRatio*b.poissonRatio)/b.youngModulusPa);
}
double hertzPatchRadiusM(double force,double radius,double modulus){
    if(!std::isfinite(force)||force<0||!std::isfinite(radius)||!(radius>0)||!std::isfinite(modulus)||!(modulus>0))throw std::invalid_argument("Hertz patch: invalid input");
    return std::cbrt(3*force*radius/(4*modulus));
}
ImpactSolution solveImpact(const ImpactInput& in){
    if(!std::isfinite(in.restitution)||in.restitution<0||in.restitution>1)throw std::invalid_argument("impact: restitution outside [0, 1]");
    const double massA=in.massAKg,massB=in.massBKg;
    const double reduced=std::isinf(massB)?massA:massA*massB/(massA+massB);
    return {(1+in.restitution)*reduced*in.normalSpeedMPerS,
            hertzContactSeconds({in.normalSpeedMPerS,massA,massB,in.radiusM,in.a.youngModulusPa,in.a.poissonRatio,in.b.youngModulusPa,in.b.poissonRatio})};
}

ContactForce::ContactForce(double sampleRateHz,std::uint32_t seed):fs_(sampleRateHz){
    if(!std::isfinite(sampleRateHz)||!(sampleRateHz>0))throw std::invalid_argument("contact force: invalid sample rate");
    // Each contact starts at its own place on the surface.
    position_=double(seed)*1e-6;
}
void ContactForce::set(const ContactForceParameters& p,const SurfaceProfile* surface) noexcept {
    p_=p;surface_=surface;
    alpha_=p.smoothingSeconds>0?1-std::exp(-1/(p.smoothingSeconds*fs_)):1;
    if(!primed_){speed_=p.speedMPerS;normal_=p.normalForceN;primed_=true;}
}
void ContactForce::generate(std::span<float> out) noexcept {
    if(p_.kind==ContactKind::Source){
        if(p_.source){p_.source->render(out);const float s=float(p_.forceScaleN);for(auto& v:out)v*=s;}
        else std::fill(out.begin(),out.end(),0.f);
        return;
    }
    const double dt=1/fs_;
    const bool roll=p_.kind==ContactKind::Roll;
    const bool stick=!roll&&p_.stickStiffnessNPerM>0&&p_.staticFriction>p_.dynamicFriction;
    const double stickMeanPole=1-std::exp(-2*Pi*5*dt);
    for(auto& sample:out){
        speed_+=alpha_*(p_.speedMPerS-speed_);normal_+=alpha_*(p_.normalForceN-normal_);
        const double step=speed_*dt;
        position_+=step;
        // The filter width only needs to change slowly; recompute on a frame grid, not a block grid.
        if((frame_++&15)==0){
            const double patch=roll&&normal_>0?2*std::cbrt(3*normal_*p_.radiusM/(4*p_.effectiveModulusPa)):0;
            width_=std::max(step,patch);
        }
        double force=0;
        if(surface_&&speed_>0&&normal_>0){
            const double height=surface_->read(position_,width_);
            force=(roll?p_.fluctuation*normal_:p_.fluctuation*p_.dynamicFriction*normal_)*height;
        }
        if(roll&&p_.outOfRound>0&&normal_>0){
            // Flat spots: a narrow cos^8 pulse per flat, zero mean (the mean of cos^8 is 35/128).
            angle_+=step/p_.radiusM;
            if(angle_>2*Pi)angle_-=2*Pi;
            const double c=std::cos(angle_*double(p_.flats)/2),c2=c*c,c4=c2*c2;
            force+=p_.outOfRound*normal_*(c4*c4-35./128.)*128./35.;
        }
        if(roll&&force<-normal_)force=-normal_; // the contact opens: a body cannot pull on the floor
        if(stick&&normal_>0){
            // A spring of tangential stiffness k dragged at the slip speed sticks until k x reaches
            // mu_s N, then slips back to (2 mu_d - mu_s) N / k: period 2 (mu_s - mu_d) N / (k v).
            stretch_+=step;
            const double hold=p_.stickStiffnessNPerM*stretch_;
            if(hold>=p_.staticFriction*normal_)stretch_=std::max(0.,(2*p_.dynamicFriction-p_.staticFriction)*normal_/p_.stickStiffnessNPerM);
            const double spring=p_.stickStiffnessNPerM*stretch_;
            stickMean_+=stickMeanPole*(spring-stickMean_);
            const double ratio=speed_/p_.stribeckSpeedMPerS;
            force+=(spring-stickMean_)/(1+ratio*ratio);
        }
        sample=float(force);
    }
}

namespace {
constexpr std::size_t CullGridFrames=64;
struct Hit {std::uint64_t frame{},sequence{};ModalExcitation excitation{};};
[[nodiscard]] bool later(const Hit& a,const Hit& b) noexcept {return a.frame!=b.frame?a.frame>b.frame:a.sequence>b.sequence;}
struct Contact {
    std::uint64_t key{};
    std::uint32_t voice{};
    ContactForce generator{48000,0};
    bool used{},sounding{},releasing{};
};
}
struct ContactAudio::Impl {
    ContactAudioSettings s;
    ModalSynth synth;
    std::vector<std::shared_ptr<const SurfaceProfile>> surfaces;
    std::vector<Contact> contacts;
    std::vector<Hit> hits;           // binary heap, capacity fixed at construction
    std::vector<float> forces;       // contacts x maximumBlockFrames
    std::vector<float> voiceScratch; // maximumVoices x maximumBlockFrames
    std::uint64_t frame{},sequence{};
    ContactAudioStats stats;
    explicit Impl(const ContactAudioSettings& settings):s(settings),synth(settings.synth){}
    [[nodiscard]] Contact* find(std::uint64_t key) noexcept {for(auto& c:contacts)if(c.used&&c.key==key)return &c;return nullptr;}
    [[nodiscard]] const Contact* find(std::uint64_t key) const noexcept {for(const auto& c:contacts)if(c.used&&c.key==key)return &c;return nullptr;}
};
ContactAudio::ContactAudio(const ContactAudioSettings& s):impl_(std::make_unique<Impl>(s)){
    if(s.maximumContacts<1||s.maximumContacts>s.synth.maximumVoices)throw std::invalid_argument("contact audio: maximum contacts must be in 1..maximumVoices");
    if(s.maximumScheduledHits<1||s.maximumScheduledHits>1u<<20)throw std::invalid_argument("contact audio: maximum scheduled hits outside 1..2^20");
    auto& m=*impl_;
    m.contacts.assign(s.maximumContacts,Contact{});
    m.hits.reserve(s.maximumScheduledHits);
    m.forces.assign(std::size_t(s.maximumContacts)*s.synth.maximumBlockFrames,0.f);
    m.voiceScratch.assign(std::size_t(s.synth.maximumVoices)*s.synth.maximumBlockFrames,0.f);
}
ContactAudio::~ContactAudio()=default;
ContactAudio::ContactAudio(ContactAudio&&) noexcept=default;
ContactAudio& ContactAudio::operator=(ContactAudio&&) noexcept=default;
ModalModelHandle ContactAudio::registerModel(std::shared_ptr<const ModalModel> model){return impl_->synth.registerModel(std::move(model));}
SurfaceHandle ContactAudio::registerSurface(std::shared_ptr<const SurfaceProfile> surface){
    if(!surface)throw std::invalid_argument("contact audio: null surface");
    impl_->surfaces.push_back(std::move(surface));
    return static_cast<SurfaceHandle>(impl_->surfaces.size()-1);
}
bool ContactAudio::valid(const ModalExcitation& e) const noexcept {
    const auto& m=*impl_;
    const auto samples=m.synth.sampleCount(e.model);
    if(e.instance==0||samples==0||!std::isfinite(e.impulseNs)||e.impulseNs<0)return false;
    if(!std::isfinite(e.contactSeconds)||e.contactSeconds<MinimumContactSeconds||e.contactSeconds>MaximumContactSeconds)return false;
    if(!std::isfinite(e.perceptualWeight)||!(e.perceptualWeight>0))return false;
    return e.sample?*e.sample<samples:std::isfinite(e.localPosition.x)&&std::isfinite(e.localPosition.y)&&std::isfinite(e.localPosition.z);
}
bool ContactAudio::valid(const ContactUpdate& u) const noexcept {
    const auto& m=*impl_;
    const auto samples=m.synth.sampleCount(u.model);
    if(u.contact==0||samples==0)return false;
    if(u.sample?*u.sample>=samples:!(std::isfinite(u.localPosition.x)&&std::isfinite(u.localPosition.y)&&std::isfinite(u.localPosition.z)))return false;
    if(!std::isfinite(u.rampSeconds)||u.rampSeconds<0||!std::isfinite(u.perceptualWeight)||!(u.perceptualWeight>0))return false;
    const auto& f=u.force;
    if(!std::isfinite(f.speedMPerS)||f.speedMPerS<0||!std::isfinite(f.normalForceN)||f.normalForceN<0||!std::isfinite(f.fluctuation)||f.fluctuation<0||!std::isfinite(f.forceScaleN))return false;
    if(!std::isfinite(f.smoothingSeconds)||f.smoothingSeconds<0)return false;
    if(f.kind==ContactKind::Source)return f.source!=nullptr;
    if(!u.surface||*u.surface>=m.surfaces.size())return false;
    if(f.kind==ContactKind::Roll)return std::isfinite(f.radiusM)&&f.radiusM>0&&std::isfinite(f.effectiveModulusPa)&&f.effectiveModulusPa>0&&f.flats>0&&std::isfinite(f.outOfRound);
    return std::isfinite(f.stribeckSpeedMPerS)&&f.stribeckSpeedMPerS>0;
}
bool ContactAudio::impact(const ModalExcitation& e,double delaySeconds) noexcept {
    auto& m=*impl_;
    if(m.hits.size()>=m.s.maximumScheduledHits||!std::isfinite(delaySeconds)||delaySeconds<0||!valid(e)){++m.stats.hitsDropped;return false;}
    m.hits.push_back({m.frame+static_cast<std::uint64_t>(std::llround(delaySeconds*m.s.synth.sampleRateHz)),m.sequence++,e});
    std::push_heap(m.hits.begin(),m.hits.end(),later);
    return true;
}
bool ContactAudio::update(const ContactUpdate& u){
    auto& m=*impl_;
    if(!valid(u))throw std::invalid_argument("contact audio: invalid update (key, model, sample, a surface for slide and roll, a source for Source, or force values)");
    const bool sourced=u.force.kind==ContactKind::Source;
    const SurfaceProfile* surface=sourced?nullptr:m.surfaces[*u.surface].get();
    Contact* c=m.find(u.contact);
    if(!c){
        for(auto& candidate:m.contacts)if(!candidate.used){c=&candidate;break;}
        if(!c){++m.stats.contactsRefused;return false;}
        c->key=u.contact;c->used=true;c->sounding=false;
        c->generator=ContactForce(m.s.synth.sampleRateHz,seedOf(u.contact));
    }
    c->releasing=false;
    c->generator.set(u.force,surface);
    const auto& f=u.force;
    const double rank=sourced?f.forceScaleN:f.kind==ContactKind::Roll?f.fluctuation*f.normalForceN:f.fluctuation*f.dynamicFriction*f.normalForceN;
    ModalDrive d;d.instance=u.contact;d.model=u.model;d.sample=u.sample;d.localPosition=u.localPosition;
    d.forceRmsN=std::max(rank,1e-12);d.rampSeconds=u.rampSeconds;d.perceptualWeight=u.perceptualWeight;d.maximumModes=u.maximumModes;d.externalForce=true;
    const auto result=m.synth.drive(d);
    c->sounding=result.outcome!=ModalExciteOutcome::RejectedQuieter&&result.outcome!=ModalExciteOutcome::Inaudible;
    c->voice=result.voice;
    return true;
}
void ContactAudio::end(std::uint64_t contact,double ramp) noexcept {
    auto& m=*impl_;
    Contact* c=m.find(contact);
    if(!c)return;
    if(!c->sounding){c->used=false;return;}
    ModalDrive d;d.instance=contact;d.model=m.synth.voice(c->voice).model;d.sample=0u;d.forceRmsN=0;d.rampSeconds=std::max(ramp,0.);
    (void)m.synth.drive(d);
    c->releasing=true;
}
void ContactAudio::stop(std::uint64_t instance) noexcept {
    auto& m=*impl_;
    if(Contact* c=m.find(instance))c->used=false;
    m.synth.stop(instance);
}
void ContactAudio::render(std::span<float> bus,std::span<float> voices,std::size_t frames){
    auto& m=*impl_;
    const std::size_t block=m.s.synth.maximumBlockFrames,slots=m.s.synth.maximumVoices;
    if(bus.size()!=frames)throw std::invalid_argument("contact audio: bus size must equal frames");
    if(!voices.empty()&&voices.size()!=slots*frames)throw std::invalid_argument("contact audio: voice buffer must be empty or maximumVoices x frames");
    std::size_t done=0;
    while(done<frames){
        while(!m.hits.empty()&&m.hits.front().frame<=m.frame){
            std::pop_heap(m.hits.begin(),m.hits.end(),later);
            const auto result=m.synth.excite(m.hits.back().excitation);
            m.hits.pop_back();
            if(result.outcome==ModalExciteOutcome::RejectedQuieter)++m.stats.impactsRefused;
            else if(result.outcome!=ModalExciteOutcome::Inaudible)++m.stats.impactsPlayed;
        }
        std::size_t chunk=std::min(frames-done,block);
        if(!m.hits.empty())chunk=std::min<std::size_t>(chunk,std::size_t(m.hits.front().frame-m.frame));
        chunk=std::min<std::size_t>(chunk,CullGridFrames-std::size_t(m.frame%CullGridFrames)); // culling must not follow the caller's block size
        for(std::size_t i=0;i<m.contacts.size();++i){
            auto& c=m.contacts[i];
            if(!c.used)continue;
            const std::span<float> force(m.forces.data()+i*block,chunk);
            c.generator.generate(force);
            if(c.sounding&&m.synth.instanceIn(c.voice)==c.key)m.synth.setDriveForce(c.voice,force.data());
        }
        const auto target=voices.empty()?std::span<float>{}:std::span<float>(m.voiceScratch.data(),slots*chunk);
        m.synth.render(bus.subspan(done,chunk),target,chunk);
        if(!voices.empty())for(std::size_t v=0;v<slots;++v)std::copy_n(m.voiceScratch.data()+v*chunk,chunk,voices.data()+v*frames+done);
        for(auto& c:m.contacts)if(c.used&&c.releasing&&(!c.sounding||m.synth.instanceIn(c.voice)!=c.key||!m.synth.driven(c.voice))){c.used=false;c.releasing=false;}
        done+=chunk;m.frame+=chunk;
    }
}
std::uint32_t ContactAudio::activeContacts() const noexcept {
    std::uint32_t n=0;for(const auto& c:impl_->contacts)n+=c.used?1u:0u;return n;
}
std::uint64_t ContactAudio::pendingHits() const noexcept {return impl_->hits.size();}
const ContactForce* ContactAudio::force(std::uint64_t contact) const noexcept {const auto* c=impl_->find(contact);return c?&c->generator:nullptr;}
const ContactAudioStats& ContactAudio::stats() const noexcept {return impl_->stats;}
ModalSynth& ContactAudio::synth() noexcept {return impl_->synth;}
const ModalSynth& ContactAudio::synth() const noexcept {return impl_->synth;}
double ContactAudio::sampleRateHz() const noexcept {return impl_->s.synth.sampleRateHz;}
}
