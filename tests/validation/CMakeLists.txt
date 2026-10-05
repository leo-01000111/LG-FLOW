find_package(GTest REQUIRED)

add_executable(lgflow_validation_tests
    test_ghia_re100.cpp
)

target_include_directories(lgflow_validation_tests PRIVATE
    ${CMAKE_SOURCE_DIR}/src
)

target_link_libraries(lgflow_validation_tests PRIVATE
    flowcore_lib
    GTest::GTest
    GTest::Main
)

# Copy runtime DLLs (gtest.dll, gtest_main.dll) next to the exe on Windows.
# Must happen before gtest_discover_tests tries to launch the binary.
if(WIN32)
    add_custom_command(TARGET lgflow_validation_tests POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_RUNTIME_DLLS:lgflow_validation_tests>
            $<TARGET_FILE_DIR:lgflow_validation_tests>
        COMMAND_EXPAND_LISTS
    )
endif()

# Use add_test (not gtest_discover_tests) so CMake does not launch the exe at
# build time — the DLL copy above only completes at the end of the build step.
# WORKING_DIRECTORY ensures cases/ and output/ resolve relative to FlowCore/.
add_test(
    NAME GhiaRe100
    COMMAND lgflow_validation_tests
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)
set_tests_properties(GhiaRe100 PROPERTIES TIMEOUT 120)
