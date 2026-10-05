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
    test_solver_pressure_corrections.cpp
    test_momentum_solver.cpp
    test_unstructured_mesh.cpp
    test_unstructured_field.cpp
    test_wall_distance.cpp
    test_spalart_allmaras.cpp
)

target_link_libraries(lgflow_tests
    PRIVATE
        flowcore_lib
        GTest::GTest
        GTest::Main
)

if(WIN32)
    add_custom_command(TARGET lgflow_tests POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_RUNTIME_DLLS:lgflow_tests>
            $<TARGET_FILE_DIR:lgflow_tests>
        COMMAND_EXPAND_LISTS)
endif()

include(GoogleTest)
gtest_discover_tests(lgflow_tests)

add_subdirectory(validation)
