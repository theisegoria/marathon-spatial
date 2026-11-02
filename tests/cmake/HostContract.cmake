# Configures tests/cmake/host, a minimal host, then checks what its ctest would run. Builds nothing.
file(REMOVE_RECURSE "${WORK}")
execute_process(COMMAND "${CMAKE_COMMAND}" -S "${CMAKE_CURRENT_LIST_DIR}/host" -B "${WORK}" -G "${GENERATOR}"
                        -DBUILD_TESTING=ON "-DLIBRARY=${LIBRARY}" RESULT_VARIABLE configured OUTPUT_VARIABLE out ERROR_VARIABLE err)
if(NOT configured EQUAL 0)
    message(FATAL_ERROR "host did not configure:\n${out}\n${err}")
endif()
execute_process(COMMAND "${CMAKE_CTEST_COMMAND}" -N WORKING_DIRECTORY "${WORK}" OUTPUT_VARIABLE listed)
if(listed MATCHES "marathon_spatial")
    message(FATAL_ERROR "the library's tests reached the host's ctest, where they are not built:\n${listed}")
endif()
if(NOT listed MATCHES "host_probe")
    message(FATAL_ERROR "the host's own test is missing:\n${listed}")
endif()
