cmake_minimum_required(VERSION 3.20)
project(LG-Flow VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ── Compiler warnings ────────────────────────────────────────────────────────
add_compile_options(
    -Wall
    -Wextra
    -Wpedantic
    $<$<CONFIG:Release>:-O3>
    $<$<CONFIG:Release>:-march=native>
    $<$<CONFIG:Debug>:-fsanitize=address,undefined>
)
add_link_options(
    $<$<CONFIG:Debug>:-fsanitize=address,undefined>
)

# ── Dependencies ─────────────────────────────────────────────────────────────
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# ── Source files ─────────────────────────────────────────────────────────────
set(FLOWCORE_SOURCES
    src/core/Mesh.cpp
    src/core/Field.cpp
    src/core/BoundaryCondition.cpp
    src/solver/NavierStokesSolver.cpp
    src/solver/PressureSolver.cpp
    src/solver/Discretization.cpp
    src/io/MeshReader.cpp
    src/io/VTKWriter.cpp
    src/utils/Logger.cpp
    src/utils/Config.cpp
)

# ── Library target (shared by solver and tests) ───────────────────────────────
add_library(flowcore_lib STATIC ${FLOWCORE_SOURCES})

target_include_directories(flowcore_lib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)

target_link_libraries(flowcore_lib PUBLIC Eigen3::Eigen)

# ── Solver executable ─────────────────────────────────────────────────────────
add_executable(lgflow src/main.cpp)

target_link_libraries(lgflow PRIVATE flowcore_lib)

# ── Tests ─────────────────────────────────────────────────────────────────────
enable_testing()
add_subdirectory(tests)
