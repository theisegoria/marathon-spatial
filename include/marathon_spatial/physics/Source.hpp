#pragma once
// The seam between physics audio and the spatial sources: a mono block source.
#include <span>

namespace mspatial::physics {
// One owner at a time, on one thread. render replaces the block's contents, does
// not allocate and leaves silence after the end.
class BlockSource {
public:
    virtual ~BlockSource()=default;
    virtual void render(std::span<float> mono) noexcept=0;
    [[nodiscard]] virtual bool finished() const noexcept {return false;}
};
}
