# LG-Flow History

This file is an append-only engineering log for work done in this repository.

## Formatting Guidelines

Use this format for every new entry:

```md
## YYYY-MM-DD HH:MM (Europe/Warsaw) - <Phase/Task Title>
- Author: <Claude Code | Codex | User>
- Status: <committed | local-uncommitted>
- Summary:
  - <what changed>
  - <what changed>
- Files changed:
  - <path>
  - <path>
- Validation:
  - <command>
  - <result>
- Risks/TODOs:
  - <remaining issue>
```

Rules:
- Append new entries to the end of this file.
- Do not edit or delete older entries.
- Always use local timezone `Europe/Warsaw`.
- Keep summaries factual and concise.
- Include exact build/test commands that were run.

---

## 2026-03-18 20:02 (Europe/Warsaw) - Initial Scaffold Created
- Author: Claude Code
- Status: committed (`5b8bd68`)
- Summary:
  - Created LG-Flow project skeleton with CMake, src/tests/cases layout, and stub interfaces.
  - Added baseline README, `.gitignore`, and placeholder tests.
- Files changed:
  - Initial repository scaffold (see commit `5b8bd68`)
- Validation:
  - No build/test evidence recorded in commit message.
- Risks/TODOs:
  - Core functionality remained stub-only by design.

## 2026-03-18 21:20 (Europe/Warsaw) - Phase 1 Hardening + Windows README Update
- Author: Claude Code + Codex
- Status: committed (`c49fbdf`)
- Summary:
  - Hardened compiler/toolchain handling in CMake for GCC/Clang/MSVC.
  - Added solver initialization-state contract (guards for `pressure()`, `velocity()`, `step()`, `run()` before `initialize()`).
  - Improved config parsing behavior documentation and error messaging.
  - Added tests for config, logger, and solver initialization contract.
  - Added Windows dependency/build instructions to README.
- Files changed:
  - `CMakeLists.txt`
  - `README.md`
  - `src/main.cpp`
  - `src/solver/Discretization.cpp`
  - `src/solver/NavierStokesSolver.cpp`
  - `src/solver/NavierStokesSolver.hpp`
  - `src/utils/Config.cpp`
  - `src/utils/Logger.cpp`
  - `tests/CMakeLists.txt`
  - `tests/test_bc.cpp`
  - `tests/test_config.cpp`
  - `tests/test_logger.cpp`
  - `tests/test_solver_init.cpp`
- Validation:
  - `C:\Program Files\CMake\bin\cmake.exe --build build --config Debug`
  - Build succeeded (`flowcore_lib`, `lgflow`, `lgflow_tests`).
  - `C:\Program Files\CMake\bin\ctest.exe --test-dir build -C Debug --output-on-failure`
  - Result: 55/55 tests passed.
- Risks/TODOs:
  - `Config::load()` does not clear previous entries before reload.
  - One solver-init test remains placeholder for invalid mesh config path.
  - `std::filesystem::exists()` pre-check in `main.cpp` is outside `try`.

## 2026-03-18 23:45 (Europe/Warsaw) - Phase 2 Mesh Geometry Core + Cleanup Fixes
- Author: Claude Code
- Status: local-uncommitted
- Summary:
  - Implemented full Mesh geometry and topology in `Mesh.cpp`: `load()` now pre-computes all cell centres, cell volumes, face areas, face normals, and face owner-neighbour connectivity following the deterministic indexing spec (x-faces first, then y-faces).
  - Boundary and interior normals assigned per spec: left/bottom boundaries use outward-pointing negative normals; right/top boundaries and interior faces use positive normals.
  - Added `std::out_of_range` guards to all Mesh accessor methods (`getCellCenter`, `getFaceArea`, `getFaceNormal`, `getNeighbors`, `getCellVolume`).
  - Cleanup fix: `Config::load()` now clears `m_entries` at entry so reload does not retain stale keys.
  - Cleanup fix: `main.cpp` existence check moved inside the `try` block and switched to the `std::error_code` overload of `std::filesystem::exists()` (no-throw).
  - Replaced all placeholder Mesh tests with full numeric assertions covering centres, volumes, face areas, normals, neighbours, and bounds checks.
  - Replaced `Initialize_ZeroNx_Throws` placeholder with a real test that writes a temp config file (`mesh.Nx=0`) and verifies `std::invalid_argument` propagates from `Mesh::load()` through `NavierStokesSolver::initialize()`.
  - Test count grew from 55 to 77; all 77 pass.
- Files changed:
  - `src/core/Mesh.cpp`
  - `src/utils/Config.cpp`
  - `src/main.cpp`
  - `tests/test_mesh.cpp`
  - `tests/test_solver_init.cpp`
- Validation:
  - `C:\Program Files\CMake\bin\cmake.exe --build build --config Debug`
  - Build succeeded (`flowcore_lib`, `lgflow`, `lgflow_tests`), zero warnings.
  - `C:\Program Files\CMake\bin\ctest.exe --test-dir build -C Debug --output-on-failure`
  - Result: 77/77 tests passed (3.60 s).
- Risks/TODOs:
  - `Field<T>` arithmetic (`norm()`, `operator+`, `operator*`) still STUB — needed for Milestone 3 residual tracking.
  - `VTKWriter::write()` and `MeshReader::read()` still STUB — needed for output and external mesh import.
  - `Discretization` operators (divergence, gradient, laplacian) still STUB — Milestone 3 scope.
  - `BoundaryCondition::applyVelocity/applyPressure` still STUB — Milestone 3 scope.
  - Lid-driven cavity validation (Ghia et al. 1982) deferred to Milestone 4.
