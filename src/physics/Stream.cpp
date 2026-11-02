#include "marathon_spatial/physics/Stream.hpp"
#include <algorithm>
#include <atomic>
#include <bit>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace mspatial::physics {
uint32_t blockSourceCallback(void* user,float* mono,uint32_t frames) noexcept {
    auto* source=static_cast<BlockSource*>(user);
    if(source->finished())return 0;
    source->render({mono,frames});
    return frames;
}
void blockCallback(void* user,uint32_t frames,uint64_t) noexcept {static_cast<PhysicsStream*>(user)->renderBlock(frames);}
namespace {
enum class Kind : std::uint8_t {Update,Impact,End,Stop};
struct Command {
    Kind kind{Kind::Stop};
    ContactUpdate update{};
    ModalExcitation hit{};
    double seconds{};
    std::uint64_t key{};
};
}
struct PhysicsStream::Impl {
    class Voice final : public BlockSource {
    public:
        Impl* owner{};std::uint32_t slot{};
        void render(std::span<float> mono) noexcept override {owner->pull(slot,mono);}
    };
    ContactAudio* audio{};
    std::vector<Command> ring;
    std::uint32_t mask{};
    alignas(64) std::atomic<std::uint64_t> head{0}; // written by the render side
    alignas(64) std::atomic<std::uint64_t> tail{0}; // written by the control side
    std::atomic<std::uint64_t> dropped{0};
    std::vector<std::atomic<std::uint64_t>> published;
    std::vector<Voice> voices;
    std::vector<float> block,bus;
    // One ring per voice, all on the render side. Sources pull what their own rate, lead and block call for, so each voice
    // reads at its own place; the synth only ever runs forward through pump.
    std::vector<float> rings;
    std::vector<std::uint32_t> readAt,level;
    std::uint32_t slots{},maximumFrames{},ringFrames{};
    bool blockDriven{};     // Set once the engine's block callback pumps; pulls then only read.
    explicit Impl(std::size_t voiceCount):published(voiceCount){}

    bool push(const Command& c) noexcept {
        const auto t=tail.load(std::memory_order_relaxed);
        if(t-head.load(std::memory_order_acquire)>mask){dropped.fetch_add(1,std::memory_order_relaxed);return false;}
        ring[t&mask]=c;
        tail.store(t+1,std::memory_order_release);
        return true;
    }
    void drain() noexcept {
        auto h=head.load(std::memory_order_relaxed);
        const auto t=tail.load(std::memory_order_acquire);
        for(;h!=t;++h){
            const auto& c=ring[h&mask];
            switch(c.kind){
            case Kind::Update:(void)audio->update(c.update);break; // validated when posted
            case Kind::Impact:(void)audio->impact(c.hit,c.seconds);break;
            case Kind::End:audio->end(c.key,c.seconds);break;
            case Kind::Stop:audio->stop(c.key);break;
            }
        }
        head.store(h,std::memory_order_release);
    }
    void publish() noexcept {
        for(std::uint32_t v=0;v<slots;++v)published[v].store(audio->synth().instanceIn(v),std::memory_order_relaxed);
    }
    void pump(std::uint32_t frames) noexcept {
        drain();
        audio->render(std::span<float>(bus.data(),frames),std::span<float>(block.data(),std::size_t(slots)*frames),frames);
        publish();
        for(std::uint32_t v=0;v<slots;++v){
            float* ring=rings.data()+std::size_t(v)*ringFrames;
            const float* from=block.data()+std::size_t(v)*frames;
            const std::uint32_t at=(readAt[v]+level[v])%ringFrames;
            const std::uint32_t first=std::min(frames,ringFrames-at);
            std::copy_n(from,first,ring+at);
            std::copy_n(from+first,frames-first,ring);
            level[v]+=frames;
            if(level[v]>ringFrames){ // Nobody is reading this voice: the oldest goes, which also keeps it from sounding late.
                const std::uint32_t lost=level[v]-ringFrames;
                readAt[v]=(readAt[v]+lost)%ringFrames;
                level[v]=ringFrames;
            }
        }
    }
    void pull(std::uint32_t slot,std::span<float> mono) noexcept {
        const auto frames=static_cast<std::uint32_t>(std::min<std::size_t>(mono.size(),ringFrames));
        if(!blockDriven)while(level[slot]<frames)pump(std::min(frames-level[slot],maximumFrames));
        const std::uint32_t copied=std::min(frames,level[slot]);
        const float* ring=rings.data()+std::size_t(slot)*ringFrames;
        const std::uint32_t first=std::min(copied,ringFrames-readAt[slot]);
        std::copy_n(ring+readAt[slot],first,mono.data());
        std::copy_n(ring,copied-first,mono.data()+first);
        readAt[slot]=(readAt[slot]+copied)%ringFrames;
        level[slot]-=copied;
        std::fill(mono.begin()+copied,mono.end(),0.f);
    }
};
PhysicsStream::PhysicsStream(ContactAudio& audio,std::uint32_t capacity){
    const auto& s=audio.synth().settings();
    if(capacity<2||capacity>1u<<20)throw std::invalid_argument("physics stream: command capacity outside 2..2^20");
    impl_=std::make_unique<Impl>(s.maximumVoices);
    auto& m=*impl_;
    m.audio=&audio;m.slots=s.maximumVoices;m.maximumFrames=s.maximumBlockFrames;
    capacity=std::bit_ceil(capacity);
    m.ring.resize(capacity);m.mask=capacity-1;
    m.voices.resize(m.slots);
    for(std::uint32_t v=0;v<m.slots;++v){m.voices[v].owner=&m;m.voices[v].slot=v;}
    m.block.assign(std::size_t(m.slots)*m.maximumFrames,0.f);
    m.bus.assign(m.maximumFrames,0.f);
    m.ringFrames=4*m.maximumFrames;
    m.rings.assign(std::size_t(m.slots)*m.ringFrames,0.f);
    m.readAt.assign(m.slots,0);
    m.level.assign(m.slots,0);
}
PhysicsStream::~PhysicsStream()=default;
bool PhysicsStream::post(const ContactUpdate& u) noexcept {
    if(!impl_->audio->valid(u)){impl_->dropped.fetch_add(1,std::memory_order_relaxed);return false;}
    Command c;c.kind=Kind::Update;c.update=u;return impl_->push(c);
}
bool PhysicsStream::postImpact(const ModalExcitation& e,double delaySeconds) noexcept {
    if(!impl_->audio->valid(e)||!std::isfinite(delaySeconds)||delaySeconds<0){impl_->dropped.fetch_add(1,std::memory_order_relaxed);return false;}
    Command c;c.kind=Kind::Impact;c.hit=e;c.seconds=delaySeconds;return impl_->push(c);
}
bool PhysicsStream::postEnd(std::uint64_t contact,double rampSeconds) noexcept {Command c;c.kind=Kind::End;c.key=contact;c.seconds=rampSeconds;return impl_->push(c);}
bool PhysicsStream::postStop(std::uint64_t instance) noexcept {Command c;c.kind=Kind::Stop;c.key=instance;return impl_->push(c);}
std::uint64_t PhysicsStream::instanceIn(std::uint32_t slot) const noexcept {return slot<impl_->slots?impl_->published[slot].load(std::memory_order_relaxed):0;}
std::uint64_t PhysicsStream::droppedCommands() const noexcept {return impl_->dropped.load(std::memory_order_relaxed);}
BlockSource& PhysicsStream::voice(std::uint32_t slot){
    if(slot>=impl_->slots)throw std::invalid_argument("physics stream: voice slot out of range");
    return impl_->voices[slot];
}
void PhysicsStream::renderBlock(std::uint32_t frames) noexcept {
    auto& m=*impl_;
    m.blockDriven=true;
    do{
        const auto n=std::min(frames,m.maximumFrames);
        m.pump(n);
        frames-=n;
    }while(frames>0);
}
void PhysicsStream::renderBus(std::span<float> mono) noexcept {
    auto& m=*impl_;
    m.drain();
    for(std::size_t at=0;at<mono.size();){
        const auto n=std::min<std::size_t>(mono.size()-at,m.maximumFrames);
        m.audio->render(mono.subspan(at,n),{},n);
        at+=n;
    }
    m.publish();
}
}
