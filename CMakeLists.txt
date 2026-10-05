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
    src/solver/MomentumSolver.cpp
    src/solver/PressureSolver.cpp
    src/solver/Discretization.cpp
    src/io/MeshReader.cpp
    src/io/VTKWriter.cpp
    src/utils/Logger.cpp
    src/utils/Config.cpp
    src/mesh/UnstructuredMesh.cpp
    src/mesh/UnstructuredBoundaryCondition.cpp
    src/io/GmshReader.cpp
    src/solver/UnstructuredDiscretization.cpp
    src/solver/UnstructuredNavierStokesSolver.cpp
    src/solver/WallDistance.cpp
    src/solver/SpallartAllmaras.cpp
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

# ── Validation executable ─────────────────────────────────────────────────────
# Run from FlowCore/: lgflow_validate_lid [config] [max_iter]
add_executable(lgflow_validate_lid
    src/validate/lid_driven_cavity_validate.cpp)
target_link_libraries(lgflow_validate_lid PRIVATE flowcore_lib)

# ── Unstructured mesh validation executable ───────────────────────────────────
add_executable(lgflow_unstruct src/validate/unstructured_channel_validate.cpp)
target_link_libraries(lgflow_unstruct PRIVATE flowcore_lib)

# ── Flat plate SA validation executable ──────────────────────────────────────
add_executable(lgflow_flat_plate src/validate/flat_plate_sa_validate.cpp)
target_link_libraries(lgflow_flat_plate PRIVATE flowcore_lib)

# ── GUI executable (SFML) ────────────────────────────────────────────────────
find_package(SFML 3 COMPONENTS Graphics Window System QUIET)
if(SFML_FOUND)
    add_executable(lgflow_ui src/ui/main_ui.cpp)
    target_include_directories(lgflow_ui PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/ui)
    # Embed the source directory so the exe can locate cases/ at runtime
    # even when launched from an unrelated working directory.
    target_compile_definitions(lgflow_ui PRIVATE
        LGFLOW_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
    target_link_libraries(lgflow_ui PRIVATE
        flowcore_lib
        SFML::Graphics
        SFML::Window
        SFML::System)
    # SFML 3 on Windows: copy DLLs next to the executable.
    # $<TARGET_RUNTIME_DLLS> only captures direct link dependencies; SFML's
    # transitive deps (freetype → brotli/bz2/libpng/zlib) must be copied
    # separately via a glob over the vcpkg bin directory.
    if(WIN32)
        add_custom_command(TARGET lgflow_ui POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
                $<TARGET_RUNTIME_DLLS:lgflow_ui>
                $<TARGET_FILE_DIR:lgflow_ui>
            COMMAND_EXPAND_LISTS)
        # Copy transitive DLL dependencies that TARGET_RUNTIME_DLLS misses.
        # (SFML Graphics → freetype → brotli, bz2, libpng, zlib)
        # A small CMake -P script avoids the file(GLOB)-vs-generator-expression
        # evaluation-order problem: the script runs at build time with the
        # correct config-specific vcpkg bin directory passed as a variable.
        if(DEFINED VCPKG_INSTALLED_DIR)
            file(WRITE "${CMAKE_BINARY_DIR}/copy_vcpkg_dlls.cmake" [[
file(GLOB dlls "${VCPKG_BIN}/*.dll")
foreach(dll ${dlls})
    get_filename_component(name "${dll}" NAME)
    set(dest "${DEST}/${name}")
    if(NOT EXISTS "${dest}" OR "${dll}" IS_NEWER_THAN "${dest}")
        file(COPY_FILE "${dll}" "${dest}" ONLY_IF_DIFFERENT)
    endif()
endforeach()
]])
            add_custom_command(TARGET lgflow_ui POST_BUILD
                COMMAND ${CMAKE_COMMAND}
                    "-DVCPKG_BIN=$<IF:$<CONFIG:Debug>,${VCPKG_INSTALLED_DIR}/x64-windows/debug/bin,${VCPKG_INSTALLED_DIR}/x64-windows/bin>"
                    "-DDEST=$<TARGET_FILE_DIR:lgflow_ui>"
                    -P "${CMAKE_BINARY_DIR}/copy_vcpkg_dlls.cmake"
                COMMENT "Copying vcpkg transitive DLLs ($<CONFIG>)")
        endif()
    endif()
    message(STATUS "LG-Flow: SFML UI target enabled (lgflow_ui)")
else()
    message(STATUS "LG-Flow: SFML not found — GUI target skipped")
endif()

# ── Tests ─────────────────────────────────────────────────────────────────────
enable_testing()
add_subdirectory(tests)
