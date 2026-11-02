// SPDX-License-Identifier: Apache-2.0
// Counts allocation and blocking calls made by one watched thread, by interposing the system functions. The test proves the
// probe works by making a deliberate allocation and lock under it first.
#include <dispatch/dispatch.h>
#include <malloc/malloc.h>
#include <os/lock.h>
#include <pthread.h>
#include <semaphore.h>
#include <unistd.h>

#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <ctime>

#include "RtProbe.h"

namespace {
std::atomic<pthread_t> gWatched{nullptr};
std::atomic<unsigned long long> gCount{0};
inline void note() {
    if (pthread_equal(pthread_self(), gWatched.load(std::memory_order_relaxed))) gCount.fetch_add(1, std::memory_order_relaxed);
}
}  // namespace

extern "C" {
// The allocator's own logging hook. It sees allocations made inside system libraries, whose calls interposing can miss.
extern void (*malloc_logger)(uint32_t, uintptr_t, uintptr_t, uintptr_t, uintptr_t, uint32_t);
static void logged(uint32_t, uintptr_t, uintptr_t, uintptr_t, uintptr_t, uint32_t) { note(); }

#ifdef MS_RT_PROBE_OFF
void rt_probe_watch_this_thread(void) {}
void rt_probe_stop(void) {}
#else
void rt_probe_watch_this_thread(void) {
    gWatched.store(pthread_self());
    malloc_logger = logged;
}
void rt_probe_stop(void) {
    malloc_logger = nullptr;
    gWatched.store(nullptr);
}
#endif
unsigned long long rt_probe_count(void) { return gCount.load(); }
void rt_probe_reset(void) { gCount.store(0); }

// Sanitizers own these same functions, so their builds compile the probe out with MS_RT_PROBE_OFF.
#ifndef MS_RT_PROBE_OFF
#define INTERPOSE(replacement, original)                                                                                         \
    __attribute__((used)) static const struct {                                                                                   \
        const void* r;                                                                                                            \
        const void* o;                                                                                                            \
    } interpose_##original __attribute__((section("__DATA,__interpose"))) = {reinterpret_cast<const void*>(replacement), reinterpret_cast<const void*>(original)}

static void* p_malloc(size_t n) { note(); return malloc(n); }
static void* p_calloc(size_t a, size_t b) { note(); return calloc(a, b); }
static void* p_realloc(void* p, size_t n) { note(); return realloc(p, n); }
static void p_free(void* p) { note(); free(p); }
static int p_posix_memalign(void** p, size_t a, size_t n) { note(); return posix_memalign(p, a, n); }
static void* p_aligned_alloc(size_t a, size_t n) { note(); return aligned_alloc(a, n); }
static int p_mutex(pthread_mutex_t* m) { note(); return pthread_mutex_lock(m); }
static void p_unfair(os_unfair_lock_t l) { note(); os_unfair_lock_lock(l); }
static int p_cond(pthread_cond_t* c, pthread_mutex_t* m) { note(); return pthread_cond_wait(c, m); }
static int p_nanosleep(const timespec* a, timespec* b) { note(); return nanosleep(a, b); }
static int p_usleep(useconds_t u) { note(); return usleep(u); }
static void p_dispatch_sync(dispatch_queue_t q, dispatch_block_t b) { note(); dispatch_sync(q, b); }
static long p_semwait(dispatch_semaphore_t s, dispatch_time_t t) { note(); return dispatch_semaphore_wait(s, t); }
INTERPOSE(p_malloc, malloc);
INTERPOSE(p_calloc, calloc);
INTERPOSE(p_realloc, realloc);
INTERPOSE(p_free, free);
INTERPOSE(p_posix_memalign, posix_memalign);
INTERPOSE(p_aligned_alloc, aligned_alloc);
INTERPOSE(p_mutex, pthread_mutex_lock);
INTERPOSE(p_unfair, os_unfair_lock_lock);
INTERPOSE(p_cond, pthread_cond_wait);
INTERPOSE(p_nanosleep, nanosleep);
INTERPOSE(p_usleep, usleep);
INTERPOSE(p_dispatch_sync, dispatch_sync);
INTERPOSE(p_semwait, dispatch_semaphore_wait);

#if __has_include(<malloc/_malloc_type.h>)
void* malloc_type_malloc(size_t, unsigned long long);
void* malloc_type_calloc(size_t, size_t, unsigned long long);
void* malloc_type_realloc(void*, size_t, unsigned long long);
int malloc_type_posix_memalign(void**, size_t, size_t, unsigned long long);
static void* p_tmalloc(size_t n, unsigned long long t) { note(); return malloc_type_malloc(n, t); }
static void* p_tcalloc(size_t a, size_t b, unsigned long long t) { note(); return malloc_type_calloc(a, b, t); }
static void* p_trealloc(void* p, size_t n, unsigned long long t) { note(); return malloc_type_realloc(p, n, t); }
static int p_tmemalign(void** p, size_t a, size_t n, unsigned long long t) { note(); return malloc_type_posix_memalign(p, a, n, t); }
INTERPOSE(p_tmalloc, malloc_type_malloc);
INTERPOSE(p_tcalloc, malloc_type_calloc);
INTERPOSE(p_trealloc, malloc_type_realloc);
INTERPOSE(p_tmemalign, malloc_type_posix_memalign);
#endif
#endif  // MS_RT_PROBE_OFF
}
