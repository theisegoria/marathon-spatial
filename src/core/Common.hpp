// SPDX-License-Identifier: Apache-2.0
// Small shared pieces: math types, aligned memory, the lock-free boundary.
#pragma once
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <new>
#include <utility>

#include "marathon_spatial/marathon_spatial.h"

namespace mspatial::core {

inline constexpr float kPi = 3.14159265358979323846f;
inline constexpr std::size_t kCacheLine = 128;  // Apple silicon's line size.

struct Vec3 {
    float x{}, y{}, z{};
    Vec3() = default;
    Vec3(float xx, float yy, float zz) : x(xx), y(yy), z(zz) {}
    Vec3(ms_vec3 v) : x(v.x), y(v.y), z(v.z) {}
    Vec3 operator+(Vec3 o) const { return {x + o.x, y + o.y, z + o.z}; }
    Vec3 operator-(Vec3 o) const { return {x - o.x, y - o.y, z - o.z}; }
    Vec3 operator*(float s) const { return {x * s, y * s, z * s}; }
};
inline float dot(Vec3 a, Vec3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline Vec3 cross(Vec3 a, Vec3 b) { return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x}; }
inline float length(Vec3 v) { return std::sqrt(dot(v, v)); }
inline Vec3 normalized(Vec3 v, Vec3 fallback = {0.0f, 0.0f, -1.0f}) {
    const float l = length(v);
    return l > 1e-9f ? v * (1.0f / l) : fallback;
}

struct Quat {
    float x{}, y{}, z{}, w{1.0f};
    Quat() = default;
    Quat(float xx, float yy, float zz, float ww) : x(xx), y(yy), z(zz), w(ww) {}
    Quat(ms_quat q) : x(q.x), y(q.y), z(q.z), w(q.w) {}
};
inline Quat conjugate(Quat q) { return {-q.x, -q.y, -q.z, q.w}; }
inline Quat normalized(Quat q) {
    const float l = std::sqrt(q.x * q.x + q.y * q.y + q.z * q.z + q.w * q.w);
    return l > 1e-9f ? Quat{q.x / l, q.y / l, q.z / l, q.w / l} : Quat{};
}
inline Quat operator*(Quat a, Quat b) {
    return {a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y, a.w * b.y - a.x * b.z + a.y * b.w + a.z * b.x,
            a.w * b.z + a.x * b.y - a.y * b.x + a.z * b.w, a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z};
}
inline Vec3 rotate(Quat q, Vec3 v) {
    const Vec3 u{q.x, q.y, q.z};
    const Vec3 t = cross(u, v) * 2.0f;
    return v + t * q.w + cross(u, t);
}

// Zero-initialised memory aligned for NEON and to cache lines. Only used off the render side.
template <class T>
class Aligned {
public:
    Aligned() = default;
    explicit Aligned(std::size_t count) { reset(count); }
    Aligned(const Aligned&) = delete;
    Aligned& operator=(const Aligned&) = delete;
    Aligned(Aligned&& o) noexcept : data_(o.data_), size_(o.size_) { o.data_ = nullptr; o.size_ = 0; }
    Aligned& operator=(Aligned&& o) noexcept {
        std::swap(data_, o.data_);
        std::swap(size_, o.size_);
        return *this;
    }
    ~Aligned() { std::free(data_); }
    void reset(std::size_t count) {
        std::free(data_);
        data_ = nullptr;
        size_ = 0;
        if (count == 0) return;
        void* p = nullptr;
        if (posix_memalign(&p, kCacheLine, count * sizeof(T)) != 0) throw std::bad_alloc();
        std::memset(p, 0, count * sizeof(T));
        data_ = static_cast<T*>(p);
        size_ = count;
    }
    T* data() { return data_; }
    const T* data() const { return data_; }
    std::size_t size() const { return size_; }
    T& operator[](std::size_t i) { return data_[i]; }
    const T& operator[](std::size_t i) const { return data_[i]; }

private:
    T* data_{};
    std::size_t size_{};
};

// Wait-free single-producer single-consumer ring of trivially copyable items.
template <class T>
class SpscRing {
public:
    void init(std::uint32_t capacity) {
        std::uint32_t c = 2;
        while (c < capacity) c <<= 1;
        items_.reset(c);
        mask_ = c - 1;
    }
    bool push(const T& item) noexcept {
        const std::uint32_t w = write_.load(std::memory_order_relaxed);
        if (w - read_.load(std::memory_order_acquire) > mask_) return false;
        items_[w & mask_] = item;
        write_.store(w + 1, std::memory_order_release);
        return true;
    }
    // All of them or none, made visible to the consumer in one store, so a consumer can never see part of a batch.
    bool pushAll(const T* items, std::uint32_t count) noexcept {
        const std::uint32_t w = write_.load(std::memory_order_relaxed);
        if (count > mask_ + 1 - (w - read_.load(std::memory_order_acquire))) return false;
        for (std::uint32_t i = 0; i < count; ++i) items_[(w + i) & mask_] = items[i];
        write_.store(w + count, std::memory_order_release);
        return true;
    }
    std::uint32_t capacity() const noexcept { return mask_ + 1; }
    bool pop(T& item) noexcept {
        const std::uint32_t r = read_.load(std::memory_order_relaxed);
        if (r == write_.load(std::memory_order_acquire)) return false;
        item = items_[r & mask_];
        read_.store(r + 1, std::memory_order_release);
        return true;
    }
    std::uint32_t pending() const noexcept {
        return write_.load(std::memory_order_acquire) - read_.load(std::memory_order_acquire);
    }

private:
    Aligned<T> items_;
    std::uint32_t mask_{};
    alignas(kCacheLine) std::atomic<std::uint32_t> write_{0};
    alignas(kCacheLine) std::atomic<std::uint32_t> read_{0};
};

// Three copies of a value: the writer fills one, the reader holds one, the third is the newest finished. Publishing and
// acquiring are one atomic exchange each, so neither side ever waits and the reader always sees a whole snapshot.
class TripleIndex {
public:
    std::uint32_t writeIndex() const noexcept { return write_; }
    std::uint32_t readIndex() const noexcept { return read_; }
    void publish() noexcept {
        const std::uint32_t old = middle_.exchange(write_ | kFresh, std::memory_order_acq_rel);
        write_ = old & kIndexMask;
    }
    bool acquire() noexcept {
        if (!(middle_.load(std::memory_order_acquire) & kFresh)) return false;
        const std::uint32_t old = middle_.exchange(read_, std::memory_order_acq_rel);
        read_ = old & kIndexMask;
        return true;
    }

private:
    static constexpr std::uint32_t kFresh = 4, kIndexMask = 3;
    alignas(kCacheLine) std::atomic<std::uint32_t> middle_{1};
    std::uint32_t write_{0};
    alignas(kCacheLine) std::uint32_t read_{2};
};

}  // namespace mspatial::core
