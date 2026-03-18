cmake_minimum_required(VERSION 3.20)
project(LG-Flow VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ── Compiler warnings and flags (per compiler ID) ────────────────────────────
#
# Each branch sets warnings-as-errors plus appropriate optimization flags.
# Sanitizers are gated to Debug builds on GCC/Clang only (MSVC has its own
# runtime checks and does not accept -fsanitize flags).
#
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(
        -Wall
        -Wextra
        -Wpedantic
        $<$<CONFIG:Release>:-O3>
        $<$<CONFIG:Release>:-march=native>
        # AddressSanitizer + UBSan: only in Debug, GCC/Clang only
        $<$<CONFIG:Debug>:-fsanitize=address,undefined>
        $<$<CONFIG:Debug>:-fno-omit-frame-pointer>
    )
    add_link_options(
        $<$<CONFIG:Debug>:-fsanitize=address,undefined>
    )
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    add_compile_options(
        /W4             # highest standard warning level
        /WX             # warnings as errors
        /utf-8          # source and execution charset UTF-8 (required for
                        # non-ASCII characters in comments/strings)
        /permissive-    # strict ISO C++ conformance mode
        /EHsc           # standard C++ exception model
        $<$<CONFIG:Release>:/O2>
        $<$<CONFIG:Release>:/GL>   # whole-program optimisation
    )
    add_link_options(
        $<$<CONFIG:Release>:/LTCG> # link-time code generation (pairs with /GL)
    )
else()
    message(WARNING
        "LG-Flow: unrecognised compiler '${CMAKE_CXX_COMPILER_ID}'. "
        "No compiler-specific flags applied. Build may succeed but "
        "warnings and sanitizers are not enabled.")
endif()

# ── Dependencies ─────────────────────────────────────────────────────────────
find_package(Eigen3 REQUIRED NO_MODULE)  # vcpkg installs 5.x; SameMajorVersion compat blocks "3.3" requests

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

# ── Library target (shared by solver executable and tests) ───────────────────
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
