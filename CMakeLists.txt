# Marathon Spatial: spatial and physics-based audio for any headphones. macOS 27+, Apple silicon M3+.
# Stands alone, or joins a host through add_subdirectory. It takes nothing from any host project.
cmake_minimum_required(VERSION 3.28)
project(MarathonSpatial VERSION 0.1.0 LANGUAGES C CXX)
set(MARATHON_SPATIAL_ROOT "${CMAKE_CURRENT_SOURCE_DIR}")

# The one switch for every test of the library. Off in a host, where they stay out of its build and out of its ctest.
option(MARATHON_SPATIAL_BUILD_TOOLS "Build the HRTF cook tool" ${PROJECT_IS_TOP_LEVEL})
option(MARATHON_SPATIAL_BUILD_TESTS "Build the library's tests" ${PROJECT_IS_TOP_LEVEL})
option(MARATHON_SPATIAL_BUILD_EXAMPLES "Build the offline examples" ${PROJECT_IS_TOP_LEVEL})

if(NOT APPLE OR NOT CMAKE_SYSTEM_PROCESSOR MATCHES "arm64")
    message(FATAL_ERROR "Marathon Spatial is macOS on Apple silicon only.")
endif()

if(PROJECT_IS_TOP_LEVEL)
    set(CMAKE_CXX_EXTENSIONS OFF)
    set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()

# Only the sources the SOFA reader needs. See third_party/libmysofa/PROVENANCE.md.
add_library(marathon_spatial_mysofa STATIC
    third_party/libmysofa/src/hdf/btree.c third_party/libmysofa/src/hdf/dataobject.c third_party/libmysofa/src/hdf/fractalhead.c
    third_party/libmysofa/src/hdf/gcol.c third_party/libmysofa/src/hdf/gunzip.c third_party/libmysofa/src/hdf/superblock.c
    third_party/libmysofa/src/hrtf/check.c third_party/libmysofa/src/hrtf/reader.c third_party/libmysofa/src/hrtf/spherical.c
    third_party/libmysofa/src/hrtf/tools.c)
target_include_directories(marathon_spatial_mysofa PUBLIC third_party/libmysofa/src/hrtf)
target_compile_options(marathon_spatial_mysofa PRIVATE -w)
target_compile_definitions(marathon_spatial_mysofa PRIVATE _GNU_SOURCE)
target_link_libraries(marathon_spatial_mysofa PUBLIC z)

# The whole engine except device output: anything linking this alone cannot open an audio device.
add_library(marathon_spatial_core STATIC
    src/core/Api.cpp src/core/Ambisonics.cpp src/core/Dsp.cpp src/core/Engine.cpp src/core/Speakers.cpp src/core/DeviceChannels.cpp src/core/Crosstalk.cpp
    src/hrtf/HrtfBuild.cpp src/hrtf/HrtfFile.cpp src/hrtf/SphereHead.cpp)
target_compile_features(marathon_spatial_core PUBLIC cxx_std_20)
target_include_directories(marathon_spatial_core PUBLIC include PRIVATE src)
target_compile_options(marathon_spatial_core PRIVATE -Wall -Wextra)
target_compile_options(marathon_spatial_core PRIVATE -fvectorize -fno-math-errno)
target_link_libraries(marathon_spatial_core PRIVATE marathon_spatial_mysofa
    "-framework Accelerate" "-framework AudioToolbox" "-framework CoreFoundation")
target_link_libraries(marathon_spatial_core INTERFACE "-framework Accelerate" "-framework AudioToolbox" "-framework CoreFoundation")

add_library(marathon_spatial_output STATIC src/core/CoreAudioOutput.cpp src/core/CameraHeadTracker.mm)
set_source_files_properties(src/core/CameraHeadTracker.mm PROPERTIES COMPILE_OPTIONS "-fobjc-arc")
target_compile_features(marathon_spatial_output PUBLIC cxx_std_20)
target_compile_options(marathon_spatial_output PRIVATE -Wall -Wextra)
target_link_libraries(marathon_spatial_output PUBLIC marathon_spatial_core "-framework CoreAudio" "-framework AVFoundation" "-framework Vision"
    "-framework CoreMedia" "-framework CoreVideo" "-framework QuartzCore")

add_library(marathon_spatial INTERFACE)
target_link_libraries(marathon_spatial INTERFACE marathon_spatial_core marathon_spatial_output)
add_library(MarathonSpatial::spatial ALIAS marathon_spatial)
add_library(MarathonSpatial::core ALIAS marathon_spatial_core)

if(MARATHON_SPATIAL_BUILD_TOOLS)
    add_executable(marathon_spatial_hrtf_cook tools/hrtf_cook.cpp)
    target_include_directories(marathon_spatial_hrtf_cook PRIVATE src)
    target_compile_features(marathon_spatial_hrtf_cook PRIVATE cxx_std_20)
    target_link_libraries(marathon_spatial_hrtf_cook PRIVATE marathon_spatial_core)
endif()

if(MARATHON_SPATIAL_BUILD_TESTS)
    enable_testing()
    # A dylib, because interposing only works from one. It counts allocations and blocking calls on one chosen thread.
    add_library(marathon_spatial_rt_probe SHARED tests/core/RtProbe.cpp)
    target_compile_features(marathon_spatial_rt_probe PRIVATE cxx_std_20)
    add_executable(marathon_spatial_core_tests tests/core/CoreTests.cpp)
    target_include_directories(marathon_spatial_core_tests PRIVATE src tests/core)
    target_compile_features(marathon_spatial_core_tests PRIVATE cxx_std_20)
    target_link_libraries(marathon_spatial_core_tests PRIVATE marathon_spatial_core marathon_spatial_rt_probe)
    target_compile_definitions(marathon_spatial_core_tests PRIVATE MS_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data/hrtf")
    add_test(NAME marathon_spatial_core COMMAND marathon_spatial_core_tests)
    set_tests_properties(marathon_spatial_core PROPERTIES LABELS "cpu;audio" TIMEOUT 600)
    # Needs a SOFA file named by MS_TEST_SOFA. None ships with the library, so without one ctest says Skipped, not Passed.
    add_test(NAME marathon_spatial_sofa COMMAND marathon_spatial_core_tests --sofa)
    set_tests_properties(marathon_spatial_sofa PROPERTIES LABELS "cpu;audio" TIMEOUT 120 SKIP_RETURN_CODE 77)
endif()

# The other parts of the library add themselves from cmake/.
file(GLOB marathon_spatial_parts CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/cmake/*.cmake")
foreach(part IN LISTS marathon_spatial_parts)
    include("${part}")
endforeach()

add_subdirectory(src/physics)
add_library(MarathonSpatial::physics ALIAS marathon_spatial_physics)
add_library(MarathonSpatial::acoustics ALIAS marathon_spatial_acoustics)
if(MARATHON_SPATIAL_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

if(MARATHON_SPATIAL_BUILD_TESTS)
    # The parts together, through the public C API, the way a host drives them.
    add_executable(marathon_spatial_integration_tests tests/integration/IntegrationTests.cpp)
    target_include_directories(marathon_spatial_integration_tests PRIVATE tests/core)
    target_link_libraries(marathon_spatial_integration_tests PRIVATE
        marathon_spatial_core marathon_spatial_physics marathon_spatial_acoustics marathon_spatial_rt_probe)
    add_test(NAME marathon_spatial_integration COMMAND marathon_spatial_integration_tests)
    set_tests_properties(marathon_spatial_integration PROPERTIES LABELS "cpu;audio" TIMEOUT 300)
    # What a host gets from add_subdirectory(... EXCLUDE_FROM_ALL), checked by configuring one.
    add_test(NAME marathon_spatial_host_contract
        COMMAND ${CMAKE_COMMAND} -DLIBRARY=${CMAKE_CURRENT_SOURCE_DIR} -DWORK=${CMAKE_CURRENT_BINARY_DIR}/host-contract
                -DGENERATOR=${CMAKE_GENERATOR} -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/cmake/HostContract.cmake)
    set_tests_properties(marathon_spatial_host_contract PROPERTIES LABELS "cpu" TIMEOUT 180)
endif()
