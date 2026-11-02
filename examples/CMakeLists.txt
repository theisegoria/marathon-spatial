# The examples render offline to WAV files. None of them links the live output library, so none can open a device.
# Add with add_subdirectory(examples) from the library's root.
add_executable(struck_plate struck_plate.cpp)
target_link_libraries(struck_plate PRIVATE marathon_spatial_physics)

if(TARGET marathon_spatial_core)
    foreach(example orbit_offline bring_your_own_audio)
        add_executable(${example} ${example}.c)
        set_target_properties(${example} PROPERTIES LINKER_LANGUAGE CXX)
        target_link_libraries(${example} PRIVATE marathon_spatial_core)
    endforeach()
    add_executable(spatial_impact spatial_impact.cpp)
    target_link_libraries(spatial_impact PRIVATE marathon_spatial_core marathon_spatial_physics)
endif()
