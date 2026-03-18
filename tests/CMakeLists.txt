find_package(GTest REQUIRED)

add_executable(lgflow_tests
    test_mesh.cpp
    test_field.cpp
    test_bc.cpp
)

target_link_libraries(lgflow_tests
    PRIVATE
        flowcore_lib
        GTest::GTest
        GTest::Main
)

include(GoogleTest)
gtest_discover_tests(lgflow_tests)
