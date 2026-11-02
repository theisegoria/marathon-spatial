# Marathon Spatial physics audio: a static library with a CPU test, built offline with no audio device.
# Add with add_subdirectory(src/physics) from the library's root, or configure this directory alone.
cmake_minimum_required(VERSION 3.28)
if(NOT PROJECT_NAME)
    project(marathon_spatial_physics LANGUAGES C CXX)
    include(CTest)
    set(MARATHON_SPATIAL_BUILD_TESTS ${BUILD_TESTING})
endif()
get_filename_component(_mspatial_root "${CMAKE_CURRENT_LIST_DIR}/../.." ABSOLUTE)
add_library(marathon_spatial_physics STATIC
    AudioAsset.cpp Contact.cpp Fracture.cpp Materials.cpp ModalAnalysis.cpp ModalAnalytic.cpp ModalSynthesis.cpp
    PhysicsCApi.cpp SampleLayers.cpp Stream.cpp)
target_include_directories(marathon_spatial_physics PUBLIC "${_mspatial_root}/include")
target_compile_features(marathon_spatial_physics PUBLIC cxx_std_23)
target_link_libraries(marathon_spatial_physics PUBLIC "-framework Accelerate" "-framework AudioToolbox" "-framework CoreFoundation")
# Eigen-analysis and resonator arithmetic are evidence: no fast-math, no contraction.
target_compile_options(marathon_spatial_physics PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-Wall -Wextra -Wpedantic -fno-fast-math -ffp-contract=off>)
if(MARATHON_SPATIAL_BUILD_TESTS)
    add_executable(marathon_spatial_physics_tests "${_mspatial_root}/tests/physics/physics_tests.cpp")
    target_link_libraries(marathon_spatial_physics_tests PRIVATE marathon_spatial_physics)
    target_compile_options(marathon_spatial_physics_tests PRIVATE -Wall -Wextra -Wpedantic -fno-fast-math -ffp-contract=off)
    add_test(NAME marathon_spatial_physics COMMAND marathon_spatial_physics_tests)
    set_tests_properties(marathon_spatial_physics PROPERTIES TIMEOUT 300 LABELS "cpu;audio")
    add_executable(marathon_spatial_physics_capi_test "${_mspatial_root}/tests/physics/physics_capi_test.c")
    set_target_properties(marathon_spatial_physics_capi_test PROPERTIES LINKER_LANGUAGE CXX)
    target_link_libraries(marathon_spatial_physics_capi_test PRIVATE marathon_spatial_physics)
    add_test(NAME marathon_spatial_physics_capi COMMAND marathon_spatial_physics_capi_test)
    set_tests_properties(marathon_spatial_physics_capi PROPERTIES TIMEOUT 60 LABELS "cpu;audio")
endif()
