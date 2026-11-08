#pragma once
// ContactAudio on a live engine. One control thread posts physics facts; the
// render thread runs the synth and hands each modal voice to its own spatial
// source. Posting and pulling never lock or allocate. The control thread may go on
// registering models and surfaces while the stream runs.
//
//   PhysicsStream stream(audio);
//   ms_engine_set_block_callback(engine,blockCallback,&stream);
//   for(slot...) ms_source_set_callback(engine,source[slot],blockSourceCallback,&stream.voice(slot));
//   game thread:  stream.post(update); stream.postImpact(hit);
//                 ms_source_set_pose(engine,source[slot],poseOf(stream.instanceIn(slot)),...);
#include "marathon_spatial/physics/Contact.hpp"
#include "marathon_spatial/physics/Source.hpp"
#include <cstdint>
#include <memory>
#include <span>

namespace mspatial::physics {
// Matches ms_source_callback. user is a BlockSource*. Returns 0 once it has finished.
uint32_t blockSourceCallback(void* user,float* mono,uint32_t frames) noexcept;
class PhysicsStream;
// Matches ms_block_callback. user is a PhysicsStream*.
void blockCallback(void* user,uint32_t frames,uint64_t frame) noexcept;

class PhysicsStream {
public:
    // audio must outlive the stream. After this, audio belongs to the render side,
    // except registerModel, registerSurface and the valid checks, which the control
    // thread may still call.
    explicit PhysicsStream(ContactAudio& audio,std::uint32_t commandCapacity=2048);
    ~PhysicsStream();
    PhysicsStream(const PhysicsStream&)=delete;
    PhysicsStream& operator=(const PhysicsStream&)=delete;
    // Control side, one thread. Each returns false for an invalid call or a full
    // queue, which is counted in droppedCommands. Commands apply at the next block.
    bool post(const ContactUpdate& update) noexcept;
    bool postImpact(const ModalExcitation& excitation,double delaySeconds=0) noexcept;
    bool postEnd(std::uint64_t contact,double rampSeconds=.03) noexcept;
    bool postStop(std::uint64_t instance) noexcept;
    bool postDrive(const ModalDrive& drive) noexcept;
    // Which instance a voice slot sounds as of the last block; 0 when idle. Any thread.
    [[nodiscard]] std::uint64_t instanceIn(std::uint32_t slot) const noexcept;
    [[nodiscard]] std::uint64_t droppedCommands() const noexcept;
    // Render side, run by the engine at the start of every block. Applies commands
    // and runs the synth, so an impact sounds even while nothing is pulling.
    void renderBlock(std::uint32_t frames) noexcept;
    // Render side: one BlockSource per voice slot. Each voice reads at its own pace
    // from up to four maximumBlockFrames of its own history, so pull sizes may differ
    // from source to source and from block to block. With renderBlock running a pull
    // never advances the synth, so a source above unity rate runs short and fills with
    // silence, and one below it falls behind: give these sources doppler 0. Without
    // renderBlock the first pull that needs frames advances the synth for every voice.
    [[nodiscard]] BlockSource& voice(std::uint32_t slot);
    void renderBus(std::span<float> mono) noexcept; // one mix instead, when per-voice sources are not used
private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};
}
