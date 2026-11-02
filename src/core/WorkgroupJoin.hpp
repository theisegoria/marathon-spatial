// SPDX-License-Identifier: Apache-2.0
// Joining the device's workgroup from its own thread. A join that does not succeed is not retried: the answer will not change
// for this thread, and a refused call on every cycle is work on the real-time path for nothing.
#pragma once
#include <cerrno>

namespace mspatial::core {

struct WorkgroupJoin {
    enum class Result { NotAttempted, Joined, AlreadyMember, Failed };
    bool attempted{};

    template <class Join>
    Result step(Join&& join) noexcept {
        if (attempted) return Result::NotAttempted;
        attempted = true;
        const int result = join();
        if (result == 0) return Result::Joined;
        return result == EALREADY ? Result::AlreadyMember : Result::Failed;
    }
    void reset() noexcept { attempted = false; }  // Only while the device thread is stopped.
};

}  // namespace mspatial::core
