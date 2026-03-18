find_package(GTest REQUIRED)

add_executable(lgflow_tests
    test_mesh.cpp
    test_field.cpp
    test_discretization.cpp
    test_bc.cpp
    test_config.cpp
    test_logger.cpp
    test_solver_init.cpp
    test_pressure_solver.cpp
    test_solver_run.cpp
    test_solver_params.cpp
    test_vtk_writer.cpp
    test_solver_convection.cpp
)

target_link_libraries(lgflow_tests
    PRIVATE
        flowcore_lib
        GTest::GTest
        GTest::Main
)

include(GoogleTest)
gtest_discover_tests(lgflow_tests)
