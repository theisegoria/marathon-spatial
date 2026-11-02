// Metal 4 hardware ray tracing for the acoustics simulator, on its own command queue. Plain C++ interface: no Metal type crosses it.
#pragma once

#include "acoustics/AcousticTrace.hpp"

#include <string>

namespace mspatial::acoustics {

struct MetalTraceConfig {
    void* device = nullptr;   // An id<MTLDevice>, borrowed. Null uses the system default device.
    std::string libraryPath;  // Path to the compiled acoustic_trace metallib.
    TraceCapacity capacity;
};

// Null with a message in `error` when this machine cannot do it: Metal 4 and Apple family 9 (M3) ray tracing are required.
// Host falls back to makeCpuTraceBackend in that case.
std::unique_ptr<TraceBackend> makeMetalTraceBackend(const MetalTraceConfig& config, std::string* error);

}  // namespace mspatial::acoustics
