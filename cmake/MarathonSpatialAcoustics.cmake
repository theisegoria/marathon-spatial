# Acoustics module of Marathon Spatial. Include from the library's own CMakeLists or from any host project.
# Depends only on the macOS SDK. Defines marathon_spatial_acoustics, its metallib, and the CPU test target.
include_guard(GLOBAL)

get_filename_component(MARATHON_SPATIAL_ROOT "${CMAKE_CURRENT_LIST_DIR}/.." ABSOLUTE)

add_library(marathon_spatial_acoustics STATIC
    ${MARATHON_SPATIAL_ROOT}/src/acoustics/AcousticScene.cpp
    ${MARATHON_SPATIAL_ROOT}/src/acoustics/AcousticSimulator.cpp
    ${MARATHON_SPATIAL_ROOT}/src/acoustics/AcousticsCApi.cpp
    ${MARATHON_SPATIAL_ROOT}/src/acoustics/CpuTraceBackend.cpp
    ${MARATHON_SPATIAL_ROOT}/src/acoustics/LateReverb.cpp
    ${MARATHON_SPATIAL_ROOT}/src/acoustics/MetalTraceBackend.mm
    ${MARATHON_SPATIAL_ROOT}/src/acoustics/RoomAcoustics.cpp)
target_compile_features(marathon_spatial_acoustics PUBLIC cxx_std_20)
target_include_directories(marathon_spatial_acoustics PUBLIC ${MARATHON_SPATIAL_ROOT}/include ${MARATHON_SPATIAL_ROOT}/src)
set_source_files_properties(${MARATHON_SPATIAL_ROOT}/src/acoustics/MetalTraceBackend.mm PROPERTIES COMPILE_OPTIONS "-fobjc-arc")
target_link_libraries(marathon_spatial_acoustics PUBLIC "-framework Foundation" "-framework Metal")

# Compiled ahead of time; pass the path to MetalTraceConfig::libraryPath.
set(MARATHON_SPATIAL_ACOUSTICS_METALLIB ${CMAKE_CURRENT_BINARY_DIR}/acoustic_trace.metallib CACHE INTERNAL "Compiled acoustic ray tracing kernels")
add_custom_command(OUTPUT ${MARATHON_SPATIAL_ACOUSTICS_METALLIB}
    COMMAND xcrun metal -std=metal4.1 -O2 -o ${MARATHON_SPATIAL_ACOUSTICS_METALLIB} ${MARATHON_SPATIAL_ROOT}/shaders/acoustic_trace.metal
    DEPENDS ${MARATHON_SPATIAL_ROOT}/shaders/acoustic_trace.metal ${MARATHON_SPATIAL_ROOT}/shaders/AcousticTraceShared.h
    VERBATIM)
add_custom_target(marathon_spatial_acoustics_metallib DEPENDS ${MARATHON_SPATIAL_ACOUSTICS_METALLIB})
# Through the library, so a host that adds this project EXCLUDE_FROM_ALL still gets the kernels by linking it.
add_dependencies(marathon_spatial_acoustics marathon_spatial_acoustics_metallib)

if(MARATHON_SPATIAL_BUILD_TESTS)
    add_executable(marathon_spatial_acoustics_tests ${MARATHON_SPATIAL_ROOT}/tests/acoustics_tests.cpp)
    target_link_libraries(marathon_spatial_acoustics_tests PRIVATE marathon_spatial_acoustics)
    add_test(NAME marathon_spatial_acoustics COMMAND marathon_spatial_acoustics_tests)
    set_tests_properties(marathon_spatial_acoustics PROPERTIES LABELS "cpu;audio" TIMEOUT 600)

    # GPU test: run by hand or from the GPU queue, under the GPU lock. Not in the default CTest set.
    add_executable(marathon_spatial_acoustics_metal_tests ${MARATHON_SPATIAL_ROOT}/tests/acoustics_metal_tests.mm)
    set_source_files_properties(${MARATHON_SPATIAL_ROOT}/tests/acoustics_metal_tests.mm PROPERTIES COMPILE_OPTIONS "-fobjc-arc")
    target_link_libraries(marathon_spatial_acoustics_metal_tests PRIVATE marathon_spatial_acoustics)
    add_dependencies(marathon_spatial_acoustics_metal_tests marathon_spatial_acoustics_metallib)
endif()
