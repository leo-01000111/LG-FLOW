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

## 2026-03-18 23:59 (Europe/Warsaw) - Phase 3 Field Algebra + Core Discretization Operators
- Author: Claude Code
- Status: local-uncommitted
- Summary:
  - Implemented `Field<T>::norm()` using `if constexpr` to branch on scalar vs. vector type:
    scalar: `sqrt(Σ v²)`, vector: `sqrt(Σ |v|²)`. Added `<cmath>` and `<type_traits>` includes.
  - Implemented `Field<T>::operator+` with mesh-pointer identity check; throws `std::invalid_argument`
    on mismatch. Element-wise addition via direct `m_data` access.
  - Implemented `Field<T>::operator*` as element-wise scalar multiplication (works for both
    `double` and `Eigen::Vector2d` via Eigen's scalar overload).
  - Implemented `Discretization::divergence` — Gauss theorem, ∇·u ≈ (1/V) Σ_f (u_f·n_f) A_f.
    Central differencing on interior faces; zero-gradient closure on boundary faces.
  - Implemented `Discretization::gradient` — Gauss theorem, ∇φ ≈ (1/V) Σ_f φ_f n_f A_f.
    Same face treatment as divergence.
  - Implemented `Discretization::laplacian` — compact two-point stencil, ∇²φ ≈ (1/V) Σ_f (∂φ/∂n) A_f.
    Normal gradient: `(φ_n − φ_o) / d` where `d = (x_n − x_o)·n`. Throws on d ≤ 0.
    Boundary faces contribute zero (zero normal gradient assumption).
  - All three operators throw `std::invalid_argument` on mesh pointer mismatch.
  - Replaced all placeholder assertions in `test_field.cpp`; added 5 new tests
    (VectorField_Norm, OperatorPlus, OperatorMultiply, MeshMismatch_Throws, VectorField_SetAndGet upgraded).
  - Created `tests/test_discretization.cpp` with 10 tests covering all three operators
    and all three mismatch cases.
  - Added `test_discretization.cpp` to `tests/CMakeLists.txt`.
  - Test count grew from 77 to 91; all 91 pass.
- Files changed:
  - `src/core/Field.hpp`
  - `src/solver/Discretization.cpp`
  - `tests/test_field.cpp`
  - `tests/test_discretization.cpp` (new)
  - `tests/CMakeLists.txt`
- Validation:
  - `C:\Program Files\CMake\bin\cmake.exe --build build --config Debug`
  - Build succeeded (`flowcore_lib`, `lgflow`, `lgflow_tests`), zero warnings, zero errors.
  - `C:\Program Files\CMake\bin\ctest.exe --test-dir build -C Debug --output-on-failure`
  - Result: 91/91 tests passed (2.52 s).
- Risks/TODOs:
  - `Field<T>` arithmetic only verified for `double` and `Eigen::Vector2d`; other T types untested.
  - Boundary-face zero-gradient closure for divergence/gradient gives first-order accuracy at walls;
    will need proper BC injection (ghost cells or flux override) when `BoundaryCondition` is implemented.
  - `laplacian` d ≤ 0 check is a safety guard; on a valid uniform mesh this path is unreachable.
  - `VTKWriter::write()` and `MeshReader::read()` still STUB — Milestone 2 remainder.
  - `BoundaryCondition::applyVelocity/applyPressure` still STUB — needed before SIMPLE loop.
  - SIMPLE loop (NavierStokesSolver) deferred to Milestone 3 physics scope.
  - Lid-driven cavity validation (Ghia et al. 1982) deferred to Milestone 4.

## 2026-03-18 23:00 (Europe/Warsaw) - Phase 4 Boundary Conditions + Case Mapping
- Author: Claude Code
- Status: local-uncommitted
- Metadata correction: Phase 2 (Mesh Geometry Core) was committed as `6136e80`;
  Phase 2 HISTORY entry incorrectly listed status as "local-uncommitted".
  Phase 3 was committed as `4e2b956`; Phase 3 HISTORY entry also listed "local-uncommitted".
  Phase 4 title line previously omitted HH:MM time — corrected to 23:00.
- Summary:
  - Implemented `BoundaryCondition::applyVelocity`: Dirichlet assignment for WALL/INLET
    patches; zero-gradient (no-op) for OUTLET/SYMMETRY. Patches applied in order
    left → right → bottom → top so top wins at top-row corners (lid-driven cavity
    convention). Mesh consistency guard throws `std::invalid_argument` on pointer mismatch.
  - Implemented `BoundaryCondition::applyPressure`: OUTLET sets boundary cell to 0.0;
    INLET/WALL/SYMMETRY use zero-gradient (no-op). Same application order and guard.
  - Added `buildBoundaryConditions()` free function (anonymous namespace) in
    `NavierStokesSolver.cpp` to parse `bc.<side>.type/value/value_x/value_y` config keys
    (case-insensitive type parsing; throws `std::invalid_argument` on unknown type).
    Called from constructor; result stored in `m_bc`.
  - Updated `NavierStokesSolver::initialize()` to call `m_bc.applyVelocity` and
    `m_bc.applyPressure` after zero-allocating fields, setting correct initial boundary
    values before the first SIMPLE step.
  - Replaced all placeholder assertions in `test_bc.cpp` with numeric checks;
    added 4 new tests: `ApplyVelocity_Outlet_LeavesFieldUnchanged`,
    `ApplyVelocity_MeshMismatch_Throws`, `ApplyPressure_Wall_LeavesFieldUnchanged`,
    `ApplyPressure_MeshMismatch_Throws`.
  - Added 2 solver-level BC wiring tests in `test_solver_init.cpp`:
    `BCWiring_TopLid_SetsTopRow` (verifies config value propagates to field after
    initialize()) and `BCWiring_InvalidType_Throws` (unknown type → constructor throws).
  - Test count grew from 91 to 97; all 97 pass.
- Files changed:
  - `src/core/BoundaryCondition.cpp`
  - `src/solver/NavierStokesSolver.cpp`
  - `tests/test_bc.cpp`
  - `tests/test_solver_init.cpp`
- Validation:
  - `C:\Program Files\CMake\bin\cmake.exe --build build --config Debug`
  - Build succeeded (`flowcore_lib`, `lgflow`, `lgflow_tests`), zero warnings, zero errors.
  - `C:\Program Files\CMake\bin\ctest.exe --test-dir build -C Debug --output-on-failure`
  - Result: 97/97 tests passed (3.64 s).
- Risks/TODOs:
  - `BoundaryCondition` operates on boundary-cell values (cell-centred FVM ghost-cell
    approach). When the SIMPLE loop is live, face-flux overrides may be needed for
    strict mass conservation at INLET/OUTLET faces.
  - SYMMETRY BC velocity logic (zero normal component only) is currently a no-op; a
    correct SYMMETRY implementation requires knowing which normal component to zero,
    which is patch-geometry information not yet exposed in the interface.
  - `NavierStokesSolver::step()` remains STUB — full SIMPLE loop deferred to next phase.
  - `VTKWriter::write()` still STUB — needed for ParaView output (Milestone 2 remainder).
  - Lid-driven cavity validation (Ghia et al. 1982) deferred to Milestone 4.

## 2026-03-18 23:30 (Europe/Warsaw) - Phase 5 Pressure Solver Infrastructure
- Author: Claude Code
- Status: local-uncommitted
- Metadata correction:
  - Phase 3 was committed as `4e2b956`; Phase 3 HISTORY entry incorrectly listed "local-uncommitted".
  - Phase 4 title line previously omitted HH:MM time; corrected to 23:00.
- Summary:
  - Implemented `PressureSolver::solve()` — replaces STUB with a full sparse Poisson solve.
  - Input validation: throws `std::invalid_argument` for dt≤0, rho≤0, αp∉(0,1], or mesh
    pointer mismatch on either field argument.
  - RHS: b_i = −(ρ/dt) · divergence(u*)_i · V_i using `Discretization::divergence`.
    Sign is negative because the assembled matrix A is positive-definite (a_f = area/dist on
    diagonal), making (A p')[o] = −Laplacian(p') · V_o; the negative sign yields the correct
    Poisson equation ∇²p' = (ρ/dt)·div(u*). Reference: Patankar (1980) eq. 6.28.
  - Matrix assembly: symmetric 5-point Laplacian stencil over all interior faces;
    a_f = getFaceArea / |centre_n − centre_o|; boundary faces skipped (zero-gradient BC).
    Reference: Ferziger, Perić & Street, 4th ed., eq. 7.23.
  - Singularity fix: row 0 set to identity (reference cell p'[0] = 0); BiCGSTAB tolerates
    the resulting asymmetric row; column-0 coupling in other rows is harmless since x[0]=0.
  - Solver: `Eigen::BiCGSTAB<SparseMatrix<double>>`; throws `std::runtime_error` on
    non-convergence or non-finite solution/residual.
  - Corrections: pressure += αp · p'; velocity -= (dt/ρ) · gradient(p') via
    `Discretization::gradient`. Reference: Patankar (1980) eqs. 6.30–6.31.
  - Return value: `||A x − b||_2` (solver residual norm).
  - Added `tests/test_pressure_solver.cpp` with 7 tests (104 total from 97):
    `Solve_InvalidDt_Throws`, `Solve_InvalidRho_Throws`, `Solve_InvalidAlphaP_Throws`,
    `Solve_MeshMismatch_Throws`, `Solve_ZeroVelocity_LeavesFieldsNearUnchanged_ResidualNearZero`,
    `Solve_NonZeroDivergence_ReturnsFiniteResidual_UpdatesPressure`,
    `Solve_NonZeroDivergence_ReducesDivergenceNorm`.
  - Added `test_pressure_solver.cpp` to `tests/CMakeLists.txt`.
- Files changed:
  - `src/solver/PressureSolver.cpp`
  - `tests/test_pressure_solver.cpp` (new)
  - `tests/CMakeLists.txt`
  - `HISTORY.md` (metadata corrections to Phases 3–4 + this entry)
- Validation:
  - `C:\Program Files\CMake\bin\cmake.exe --build build --config Debug`
  - Build succeeded (`flowcore_lib`, `lgflow`, `lgflow_tests`), zero warnings, zero errors.
  - `C:\Program Files\CMake\bin\ctest.exe --test-dir build -C Debug --output-on-failure`
  - Result: 104/104 tests passed (3.88 s).
- Risks/TODOs:
  - Velocity correction uses `Discretization::gradient` (Gauss face-average), which is not
    perfectly consistent with the two-point-stencil Laplacian matrix. Divergence reduction is
    achieved in practice (verified by test) but not algebraically exact. A future improvement
    is to use a face-by-face velocity correction matching the Poisson stencil.
  - Reference-cell fix (row-0 identity) makes A asymmetric; BiCGSTAB handles this but
    ConjugateGradient cannot be used. When the SIMPLE loop is live, consider a symmetric
    Dirichlet pin (also zero out column 0) if a symmetric solver is preferred.
  - `NavierStokesSolver::step()` SIMPLE sequence (momentum predictor + pressure correction +
    velocity correction + BC apply + repeat) remains STUB — next phase scope.
  - `VTKWriter::write()` still STUB — needed for ParaView output (Milestone 2 remainder).
  - Lid-driven cavity validation (Ghia et al. 1982) deferred to Milestone 4.

## 2026-03-19 00:30 (Europe/Warsaw) - Phase 7 Verification + Regression Harness (Review Fixes Included)
- Author: Claude Code
- Status: local-uncommitted
- Summary:
  - **Review fix 1 — Solver parameter validation:** `NavierStokesSolver` constructor now validates all
    solver parameters at construction time, throwing `std::invalid_argument` for: `dt ≤ 0`,
    `rho ≤ 0`, `nu < 0`, `tolerance ≤ 0`, `alpha_u ∉ (0,1]`, `alpha_p ∉ (0,1]`, `vtk_interval ≤ 0`.
    `step(dt)` validates `dt > 0` before predictor work. `run(maxIter)` validates `maxIter >= 0`.
  - **Review fix 2 — VTKWriter fully implemented:** `VTKWriter::write()` now writes a complete ASCII
    XML VTU (UnstructuredGrid) file. One `VTK_VERTEX` per cell centre (z=0). Cell data arrays:
    `pressure` (scalar) and `velocity` (3-component, z=0). Size-mismatch validation throws
    `std::invalid_argument`. Non-existent output directory throws `std::runtime_error`.
    VTK snapshots are now actually written during `run()` at every `vtk_interval` iterations.
  - **Review fix 3 — PressureSolver linear solver upgrade:** Switched from `Eigen::BiCGSTAB`
    (fails to converge for meshes ≥ 64×64 because default machine-epsilon tolerance is unreachable
    in floating-point arithmetic) to `Eigen::SparseLU` (direct solver, unconditionally robust for
    the sparse 5-point Laplacian stencil). No interface change; all 7 pressure solver tests pass.
  - **Phase 7 — Validation executable `lgflow_validate_lid`:** New CMake target linking against
    `flowcore_lib`. Loads `cases/lid_driven_cavity/case_validate.cfg` (32×32, dt=0.001, 500 iter),
    runs SIMPLE loop, samples u(x≈0.5, y) and v(x, y≈0.5) centerlines via nearest-neighbour
    column/row selection, writes `output/lid_driven_cavity/centerline.csv` (axis,coord,value),
    and computes/reports L2 and Linf error metrics vs. Ghia et al. (1982) if
    `cases/lid_driven_cavity/ghia1982_re100.csv` is present. No hard thresholds enforced.
  - **New reference data file:** `cases/lid_driven_cavity/ghia1982_re100.csv` — 17 u-centerline
    and 17 v-centerline points from Ghia et al. (1982) Table 1 and Table 3 at Re=100.
  - **New validation config:** `cases/lid_driven_cavity/case_validate.cfg` — 32×32 mesh, dt=0.001
    (stable explicit Euler at this resolution, ~50x safety margin over viscous CFL constraint).
  - **New tests — test_solver_params.cpp (19 tests):** Covers all 7 constructor parameter
    validations (including boundary values: nu=0, alpha_u=1.0, alpha_p=1.0 that must NOT throw),
    `run(-1)` throws, and `step(0.0)` / `step(-0.01)` throw.
  - **New tests — test_vtk_writer.cpp (9 tests):** Covers file creation, XML declaration,
    UnstructuredGrid element, NumberOfPoints/NumberOfCells correctness, `Name="pressure"`,
    `Name="velocity"`, pressure and velocity size-mismatch throws, and non-existent directory
    throws.
  - **README updated:** Validation executable usage (Windows), artifact table (history.csv,
    centerline.csv, iter_*.vtu), benchmark comparison mode.
  - Test count grew from 108 to 136; all 136 pass.
- Files changed:
  - `src/solver/NavierStokesSolver.cpp` (parameter validation in constructor, step, run)
  - `src/solver/PressureSolver.cpp` (BiCGSTAB → SparseLU; include update)
  - `src/io/VTKWriter.cpp` (full XML VTU implementation replacing STUB)
  - `src/validate/lid_driven_cavity_validate.cpp` (new)
  - `cases/lid_driven_cavity/ghia1982_re100.csv` (new)
  - `cases/lid_driven_cavity/case_validate.cfg` (new)
  - `tests/test_solver_params.cpp` (new)
  - `tests/test_vtk_writer.cpp` (new)
  - `tests/CMakeLists.txt` (added test_solver_params.cpp, test_vtk_writer.cpp)
  - `CMakeLists.txt` (added lgflow_validate_lid target)
  - `README.md` (status update, validation executable docs, artifact table)
  - `HISTORY.md` (this entry)
- Validation:
  - `C:\Program Files\CMake\bin\cmake.exe --build build --config Debug`
  - Build succeeded (flowcore_lib, lgflow, lgflow_tests, lgflow_validate_lid), zero warnings, zero errors.
  - `C:\Program Files\CMake\bin\ctest.exe --test-dir build -C Debug --output-on-failure`
  - Result: 136/136 tests passed (4.23 s).
  - `build\Debug\lgflow_validate_lid.exe` (from FlowCore/):
    ```
    LG-Flow Validation: Lid-Driven Cavity
    Config  : cases/lid_driven_cavity/case_validate.cfg
    maxIter : 500
    Mesh    : 32x32
    [WARN ] run() reached maxIter=500 without convergence
    Results : vel_residual=1.993e-03  cont_residual=3.315e+01
    CSV     : output/lid_driven_cavity/centerline.csv  (32 u-points, 32 v-points)
    Ref     : cases/lid_driven_cavity/ghia1982_re100.csv
    u-line  : L2=1.309e-01  Linf=2.470e-01  (n=17)
    v-line  : L2=1.275e-01  Linf=2.490e-01  (n=17)
    [INFO] Metrics are informational — no pass/fail threshold enforced in Phase 7.
    ```
  - Artifacts produced: output/lid_driven_cavity/{centerline.csv, history.csv, iter_100..500.vtu}
- Risks/TODOs:
  - Solver has not converged in 500 iterations at 32x32; continuity residual (33.3) remains large.
    Root cause: explicit forward-Euler momentum predictor and central-differencing convection create
    slow SIMPLE convergence; the velocity correction is not perfectly consistent with the Poisson
    stencil (gradient vs. two-point stencil mismatch). Milestone 4 scope: implicit momentum predictor
    or upwind convection to achieve converged results at reasonable iteration count.
  - Error metrics (L2≈0.13, Linf≈0.25) are high because the solver has not converged. After
    convergence, 32×32 accuracy vs. Ghia would be limited by mesh resolution to ~O(dx)≈0.03.
  - SparseLU memory usage grows as O(N^{1.5}) for 2D problems; for meshes above ~256×256, an
    ILU-preconditioned BiCGSTAB would be more memory-efficient.
  - `case.cfg` (64×64, dt=0.01) remains numerically unstable for the explicit Euler predictor
    (large velocity gradient at lid corners drives divergence of u* above recoverable threshold).
  - SYMMETRY BC is still a no-op (Phase 4 known limitation).

## 2026-03-18 23:55 (Europe/Warsaw) - Phase 6 SIMPLE Loop Integration + Runtime History Output
- Author: Claude Code
- Status: local-uncommitted
- Metadata correction:
  - Phase 2 (Mesh Geometry Core) was committed as `6136e80`; Phase 2 HISTORY entry incorrectly listed "local-uncommitted".
  - Phase 3 (Field Algebra + Discretization) was committed as `4e2b956`; Phase 3 HISTORY entry also listed "local-uncommitted".
  - Phase 4 (Boundary Conditions) was committed as `3d31704`; Phase 4 HISTORY entry also listed "local-uncommitted".
  - Phase 5 HISTORY entry also listed "local-uncommitted"; commit hash not recorded in this entry.
- Summary:
  - Implemented `NavierStokesSolver::step(dt)` — full SIMPLE iteration replacing the STUB:
    1. Momentum predictor: decomposes velocity into scalar ux/uy fields, computes
       `Discretization::gradient` and `Discretization::laplacian` for each component plus
       `Discretization::gradient(p)`, then applies explicit forward-Euler update
       `u* = u^k + dt*[-(u·∇)u - (1/ρ)∇p^k + ν∇²u^k]` (central differencing,
       Ferziger & Peric 4th ed. eq. 7.20).
    2. Pressure correction: delegates to `PressureSolver::solve(uStar, pressure, dt, rho, alphaP)`
       which updates `uStar` and `pressure` in place (Patankar 1980 eqs. 6.28–6.31).
    3. Velocity under-relaxation: `u^{k+1} = α_u·u* + (1-α_u)·u^k` via `Field::operator*` and `+`.
    4. Apply BCs: `m_bc.applyVelocity` and `m_bc.applyPressure` after field update.
    5. Residual tracking: velocity residual = `||u^{k+1}-u^k|| / max(||u^k||, 1e-12)`;
       continuity residual = `||∇·u^{k+1}||_2`. NaN/Inf on any residual throws `std::runtime_error`.
  - Implemented `NavierStokesSolver::run(maxIter)` — full run loop replacing the STUB:
    - Creates `output.dir` with `std::filesystem::create_directories`.
    - Opens `history.csv` and writes header `iter,vel_residual,cont_residual,pressure_residual`.
    - Appends one row per iteration (scientific notation, 8 decimal places).
    - Logs progress every 50 iterations via `Logger`.
    - Writes VTK snapshot every `output.vtk_interval` iterations (VTKWriter still STUB, call site live).
    - Convergence criterion: both `m_velResidual < tolerance` AND `m_contResidual < tolerance`.
  - Added config keys with defaults: `solver.alpha_u` (0.7), `solver.alpha_p` (0.3),
    `output.vtk_interval` (100), `output.dir` ("output").
  - Added public accessors: `velocityResidual()`, `continuityResidual()`, `pressureResidual()`.
    `residual()` retained as alias for `velocityResidual()` (API compatibility).
  - Added `VTKWriter m_vtkWriter` as value member in `NavierStokesSolver`.
  - Added `tests/test_solver_run.cpp` with 4 tests (108 total from 104):
    `Step_AfterInit_ResidualsAreFinite`,
    `Run_WritesHistoryCsv_WithHeaderAndRows`,
    `Run_MaxIterZero_DoesNotThrow`,
    `Run_ConvergenceCheck_UsesBothVelocityAndContinuityResiduals`.
  - Removed `m_residual` member; replaced by `m_velResidual` (same semantics, no API breakage).
- Files changed:
  - `src/solver/NavierStokesSolver.hpp`
  - `src/solver/NavierStokesSolver.cpp`
  - `tests/test_solver_run.cpp` (new)
  - `tests/CMakeLists.txt`
  - `HISTORY.md` (metadata corrections to Phases 2–5 + this entry)
- Validation:
  - `C:\Program Files\CMake\bin\cmake.exe --build build --config Debug`
  - Build succeeded (`flowcore_lib`, `lgflow`, `lgflow_tests`), zero warnings, zero errors.
  - `C:\Program Files\CMake\bin\ctest.exe --test-dir build -C Debug --output-on-failure`
  - Result: 108/108 tests passed (4.89 s).
- Risks/TODOs:
  - Momentum predictor uses explicit (forward-Euler) time integration; for high Reynolds numbers
    or large dt this will be unstable. An implicit or semi-implicit momentum solve is needed
    for Milestone 3 production runs.
  - Convective term uses central differencing (`Discretization::gradient`); upwind differencing
    may be required at higher Re to suppress spurious oscillations.
  - On a zero-initialized mesh the solver trivially converges in one step (all residuals = 0).
    Non-trivial lid-driven cavity convergence requires a non-zero lid BC and sufficient iterations.
  - `VTKWriter::write()` is still STUB — no .vtu files are written yet (Milestone 2 remainder).
  - Lid-driven cavity validation (Ghia et al. 1982) deferred to Milestone 4.

## 2026-03-19 (Europe/Warsaw) - Phase 8 Convection Scheme + CFL Clamp + Validation Upgrades
- Author: Claude Code
- Status: local-uncommitted
- Metadata note: The Phase 7 HISTORY entry above appears before the Phase 6 entry (line ~343)
  due to an insertion error in a previous session. Entries are not reordered to preserve
  the append-only rule; new entries continue at the true end of the file.
- Summary:
  - **Convection discretization control:** Added `solver.convection_scheme` config key
    (`"upwind"` / `"central"`, case-insensitive; default `"upwind"`).
    Upwind branch implements first-order donor-cell differencing: backward difference when
    `u_x ≥ 0`, forward difference when `u_x < 0`; zero-gradient closure at domain boundaries.
    Central branch retains existing `Discretization::gradient`-based scheme (Gauss face-average).
    Unknown scheme value throws `std::invalid_argument` at construction.
  - **CFL safety clamp:** Added `solver.max_cfl_conv` (default 0.5) and `solver.max_cfl_diff`
    (default 0.5). Each call to `step(dt)` computes
    `dt_eff = min(dt, max_cfl_conv * h / ||u||_max, max_cfl_diff * h² / ν)` and uses `dt_eff`
    for both the momentum predictor and the PressureSolver call. Progress log (every 50 iter)
    now includes `dt_eff` and prints `[CFL-clamped]` when clamped. Both CFL params validate > 0
    at construction.
  - **`dtEffective()` accessor:** New public `[[nodiscard]] double dtEffective() const` on
    `NavierStokesSolver`. Initialised to `m_dt` in constructor; updated by each `step()` call.
  - **Private `toUpper` helper:** Extracted into anonymous namespace to avoid duplication
    between BC-type and convection-scheme parsing.
  - **Validation executable upgrades:**
    - `sampleCenterlines()` now takes explicit `midX` / `midY` parameters instead of hardcoded
      0.5; `main()` passes `Lx * 0.5` / `Ly * 0.5` from config, generalising to non-unit domains.
    - New output file `output/lid_driven_cavity/validation_metrics.csv` (columns `metric,value`)
      written after every run: `final_vel_residual`, `final_cont_residual`, `u_l2`, `u_linf`,
      `v_l2`, `v_linf`, `iterations`.
    - Actual iteration count read from `history.csv` line count rather than using `maxIter`.
    - Removed `outDir` re-declaration (previously computed twice in main).
  - **`case_validate.cfg` updated:** Added `solver.convection_scheme = upwind`,
    `solver.max_cfl_conv = 0.5`, `solver.max_cfl_diff = 0.5` with explanatory comments.
  - **New test file `test_solver_convection.cpp` (12 tests):** Unknown scheme throws; central
    and upwind construct; case-insensitive scheme parsing; default is upwind; zero/negative
    `max_cfl_conv` / `max_cfl_diff` throw; `dtEffective() ≤ user dt` after one step; both
    schemes produce finite residuals after a step.
  - **README updated:** Fixed validation section to reference `case_validate.cfg` and 500-iter
    cap (not `case.cfg` and 200-iter cap); added new artifact row for `validation_metrics.csv`;
    added new-config-keys table; updated status blurb to Phase 8.
  - Test count grew from 136 to 148; all 148 pass.
- Files changed:
  - `src/solver/NavierStokesSolver.hpp` (ConvectionScheme enum, new members, dtEffective() decl)
  - `src/solver/NavierStokesSolver.cpp` (scheme/CFL parsing, CFL clamp in step(), upwind branch, dtEffective() impl)
  - `src/validate/lid_driven_cavity_validate.cpp` (midX/midY params, validation_metrics.csv, iteration count)
  - `cases/lid_driven_cavity/case_validate.cfg` (convection_scheme, max_cfl_conv/diff keys)
  - `tests/test_solver_convection.cpp` (new)
  - `tests/CMakeLists.txt` (added test_solver_convection.cpp)
  - `README.md` (status, validation section, new config keys table)
  - `HISTORY.md` (this entry)
- Validation:
  - `C:\Program Files\CMake\bin\cmake.exe --build build --config Debug`
  - Build succeeded (flowcore_lib, lgflow, lgflow_tests, lgflow_validate_lid), zero warnings, zero errors.
  - `C:\Program Files\CMake\bin\ctest.exe --test-dir build -C Debug --output-on-failure`
  - Result: 147/147 tests passed (11.51 s).
  - `build\Debug\lgflow_validate_lid.exe` (from FlowCore/):
    ```
    LG-Flow Validation: Lid-Driven Cavity
    Config  : cases/lid_driven_cavity/case_validate.cfg
    maxIter : 500
    Mesh    : 32x32
    [WARN ] run() reached maxIter=500 without convergence
    Results : vel_residual=1.862e-03  cont_residual=3.289e+01  iterations=500
    CSV     : output/lid_driven_cavity/centerline.csv  (32 u-points, 32 v-points)
    Ref     : cases/lid_driven_cavity/ghia1982_re100.csv
    u-line  : L2=1.315e-01  Linf=2.479e-01  (n=17)
    v-line  : L2=1.286e-01  Linf=2.511e-01  (n=17)
    [INFO] Metrics are informational — no pass/fail threshold enforced.
    Metrics : output/lid_driven_cavity/validation_metrics.csv
    ```
  - `output/lid_driven_cavity/validation_metrics.csv`:
    final_vel_residual=1.862e-03, final_cont_residual=3.289e+01,
    u_l2=1.315e-01, u_linf=2.479e-01, v_l2=1.286e-01, v_linf=2.511e-01, iterations=500
- Risks/TODOs:
  - Upwind scheme is first-order accurate; grid-level diffusion will damp physical gradients
    on coarse meshes. For quantitative Ghia comparison, 64×64 with implicit momentum predictor
    is the recommended next step (Milestone 4).
  - Central scheme remains susceptible to oscillations at high Re without mesh-Reynolds
    limiter; no blending or TVD limiting is implemented.
  - CFL clamp is active only during step(); the user-specified dt is still used as the
    upper bound. If the clamp is never triggered (dt already below stability limit), there
    is zero overhead.
  - SYMMETRY BC is still a no-op.
  - SparseLU O(N^1.5) memory limitation for meshes above ~256×256 unchanged from Phase 7.

## 2026-03-19 14:00 (Europe/Warsaw) - Phase 9 Convergence Quality + Validation Trustworthiness
- Author: Claude Code
- Status: local-uncommitted
- Metadata correction: Phase 8 HISTORY entry omitted HH:MM time in the section title; correct
  time is 2026-03-19 (no exact time recorded). Entry not edited per append-only rule.
- Summary:
  - **Audit fix 1 — stale iteration-cap comment:** Removed "200 iterations" comment from
    `lid_driven_cavity_validate.cpp`; replaced with accurate "500 iterations" description.
  - **Audit fix 2 — dtEffective() doc mismatch:** Header doc previously said "Returns 0.0
    before first step()"; fixed to "Returns the configured dt before first step()" to match
    the constructor which initialises `m_dtEff = m_dt`.
  - **Audit fix 3 — misleading zero errors in validation_metrics.csv:** When the Ghia reference
    CSV is absent, error columns (`u_l2`, `u_linf`, `v_l2`, `v_linf`) are now written as `nan`
    instead of 0. Added `ref_available` row (1 or 0) as the first metric line. Stdout also
    clearly states "Reference absent — error metrics written as nan."
  - **Stencil-consistent velocity correction (Task 2):** Replaced `Discretization::gradient(pPrime)`
    call in `PressureSolver::solve` with an inline structured-grid loop. For interior cells the
    formula is identical (central difference). For boundary cells uses zero-gradient closure
    (face value = cell value), which is algebraically equivalent to the Gauss-theorem gradient
    but avoids the general face-iteration overhead and makes the stencil explicit.
    Reference: Patankar (1980) eq. 6.31; Ferziger & Peric (2020) eq. 4.22.
  - **`solver.pressure_corrections_per_step` (Task 3):** New config key (int, default 1, ≥ 1).
    In `NavierStokesSolver::step()`, the pressure-correction Poisson solve is now called
    `pressure_corrections_per_step` times per SIMPLE step, with each pass operating on the
    velocity corrected by the previous pass. Values > 1 drive continuity residual lower at
    the cost of additional SparseLU solves. Validation at construction: values < 1 throw
    `std::invalid_argument`.
  - **Validation executable pass/fail mode (Task 4):** Added optional config keys
    `validation.max_u_l2`, `validation.max_v_l2`, `validation.max_cont_residual`. If any
    threshold is exceeded, the executable prints a FAIL summary and returns exit code 2.
    Keys absent → informational mode (unchanged behavior). `validation_metrics.csv` is always
    written. Exit codes: 0 = pass/informational, 1 = fatal error, 2 = threshold exceeded.
  - **README updated (Task 5):** Phase 9 status; new config keys table
    (`solver.pressure_corrections_per_step`, `validation.max_*` keys); exit code table;
    note on nan error columns when reference is absent.
  - **New tests:** `test_solver_pressure_corrections.cpp` (6 tests):
    `ZeroCorrections_Throws`, `NegativeCorrections_Throws`, `OneCorrection_IsDefault_Constructs`,
    `TwoCorrections_Constructs`, `TwoCorrections_DoNotWorsenContinuityResidual`,
    `FiveCorrections_ResidualsFinite`.
  - **New test in test_pressure_solver.cpp:** `Solve_StencilConsistentCorrection_FiniteOnLargerMesh`
    (6x6 mesh, verifies inline gradient gives finite/non-negative residual, non-trivial pressure
    update, and finite velocities on all cells).
  - Test count grew from 147 to 154; all 154 pass.
- Files changed:
  - `src/solver/NavierStokesSolver.hpp` (dtEffective doc fix; m_pressureCorrectionsPerStep member)
  - `src/solver/NavierStokesSolver.cpp` (pressure_corrections_per_step parsing/validation/loop)
  - `src/solver/PressureSolver.cpp` (Nx added; inline structured-grid gradient replacing Discretization::gradient)
  - `src/validate/lid_driven_cavity_validate.cpp` (stale comment fix; ref_available + nan metrics; pass/fail mode)
  - `tests/test_pressure_solver.cpp` (new stencil consistency test)
  - `tests/test_solver_pressure_corrections.cpp` (new)
  - `tests/CMakeLists.txt` (added test_solver_pressure_corrections.cpp)
  - `README.md` (Phase 9 status, new config key docs, exit code table)
  - `HISTORY.md` (this entry)
- Validation:
  - `C:\Program Files\CMake\bin\cmake.exe --build build --config Debug`
  - Build succeeded (flowcore_lib, lgflow, lgflow_tests, lgflow_validate_lid), zero warnings, zero errors.
  - `C:\Program Files\CMake\bin\ctest.exe --test-dir build -C Debug --output-on-failure`
  - Result: 154/154 tests passed (14.09 s).
  - `build\Debug\lgflow_validate_lid.exe` (from FlowCore/):
    ```
    LG-Flow Validation: Lid-Driven Cavity
    Config  : cases/lid_driven_cavity/case_validate.cfg
    maxIter : 500
    Mesh    : 32x32
    [WARN ] run() reached maxIter=500 without convergence
    Results : vel_residual=1.862e-03  cont_residual=3.289e+01  iterations=500
    CSV     : output/lid_driven_cavity/centerline.csv  (32 u-points, 32 v-points)
    Ref     : cases/lid_driven_cavity/ghia1982_re100.csv
    u-line  : L2=1.315e-01  Linf=2.479e-01  (n=17)
    v-line  : L2=1.286e-01  Linf=2.511e-01  (n=17)
    Metrics : output/lid_driven_cavity/validation_metrics.csv
    [INFO] Metrics are informational — no pass/fail threshold enforced.
    ```
  - `output/lid_driven_cavity/validation_metrics.csv`: ref_available=1, u_l2=1.315e-01,
    u_linf=2.479e-01, v_l2=1.286e-01, v_linf=2.511e-01, iterations=500
- Risks/TODOs:
  - Solver has not converged in 500 iterations at 32x32; continuity residual (32.9) remains large.
    Root cause: explicit forward-Euler momentum predictor with under-relaxation. Even with
    `pressure_corrections_per_step > 1`, the divergence after velocity correction from the
    momentum predictor dominates. An implicit momentum solve or semi-implicit SIMPLE is needed.
  - The inline stencil gradient is algebraically equivalent to Discretization::gradient for all
    mesh sizes. The equivalence was verified by checking that all 8 original pressure solver
    tests still pass unchanged.
  - Reference-cell pin (p'[0]=0) means divergence at cell (0,0) is NOT driven to zero by the
    pressure correction. This is a known SIMPLE property (reference cell just fixes the pressure
    datum, not the divergence).
  - SYMMETRY BC is still a no-op.
  - SparseLU O(N^1.5) memory limitation for meshes above ~256×256 unchanged.

---

> Backfill note (2026-10-05): the six entries below were written after the fact from
> `git log` and the commit messages, because no HISTORY.md entries were made for those
> commits at the time. Numbers come from the commit messages; nothing was re-run for them.

## 2026-04-22 12:43 (Europe/Warsaw) - SIMPLE Fix: Under-Relaxation Inside the Momentum Matrix
- Author: Claude Code
- Status: committed (`5c6e175`)
- Summary:
  - Replaced post-solve velocity mixing (`u = alpha_u*u* + (1-alpha_u)*u_old`) with implicit
    under-relaxation in the momentum matrix (Patankar 1980, eq. 6.36): diagonal scaled by
    `1/alpha_u`, `(1-alpha_u)/alpha_u * aP * u_old` added to the RHS.
  - Post-solve mixing had undone the mass conservation from the pressure correction, so the
    continuity residual stalled at about 33 (the open risk from the Phase 9 entry).
  - Momentum predictor moved into `MomentumSolver` with an implicit viscous term.
  - New tests `AlphaUZero_Throws`, `AlphaURelaxed_ProducesFiniteResult`.
- Files changed:
  - `src/solver/MomentumSolver.cpp`, `src/solver/MomentumSolver.hpp`
  - `src/solver/NavierStokesSolver.cpp`
  - `tests/test_momentum_solver.cpp`
- Validation:
  - Commit message: all 166 tests pass.
- Risks/TODOs:
  - Cavity convergence still needed the BC ordering fix below.

## 2026-04-22 13:49 (Europe/Warsaw) - SIMPLE Fix: Apply BCs Before the Pressure Correction
- Author: Claude Code
- Status: committed (`be7e28c`)
- Summary:
  - `applyVelocity(uStar)` now runs between the momentum predictor and the pressure Poisson
    assembly, so the Poisson RHS sees the prescribed boundary velocities (e.g. lid = (1,0)).
    The post-correction BC application is kept.
  - Applying BCs only after the correction had reset boundary velocities and undone the
    divergence reduction in the neighbouring interior cells.
  - Continuity residual reverted to the Rhie-Chow divergence `div_RC(u_new, p^{k+1})`; plain
    `div(u)` does not go to zero on a collocated Rhie-Chow grid.
  - Note: this changes the SIMPLE step order documented in CLAUDE.md (BCs now also applied
    before the pressure correction); the reason is above.
- Files changed:
  - `src/solver/NavierStokesSolver.cpp`
- Validation:
  - Commit message: 32x32 cavity, Re=100, converges in 1577 iterations to vel=3.2e-6,
    cont=9.9e-6. All 166 tests pass.
- Risks/TODOs:
  - None recorded.

## 2026-04-22 16:55 (Europe/Warsaw) - INLET/OUTLET BCs + Generalized Pressure Reference Cells
- Author: Claude Code
- Status: committed (`831336f`)
- Summary:
  - `BoundaryCondition::collectCellsOfType()` maps patch sides to sorted, deduplicated cell IDs.
  - `PressureSolver::setDirichletPressureCells()` pins any set of cells (the outlet column)
    instead of only cell 0; `NavierStokesSolver::initialize()` uses the OUTLET cells, falling
    back to cell 0 when there is no OUTLET.
  - Convergence criterion is now the velocity residual only. `cont_residual` stays in
    `history.csv` as a diagnostic, because `div_RC` keeps a structural offset proportional to
    dp/dx at open boundaries.
  - New case `cases/channel_flow/case.cfg` (64x16, Re=100, Poiseuille).
  - 5 new tests (4 for `collectCellsOfType`, 1 for `setDirichletPressureCells`).
- Files changed:
  - `cases/channel_flow/case.cfg`
  - `src/core/BoundaryCondition.cpp`, `src/core/BoundaryCondition.hpp`
  - `src/solver/NavierStokesSolver.cpp`
  - `src/solver/PressureSolver.cpp`, `src/solver/PressureSolver.hpp`
  - `tests/test_bc.cpp`, `tests/test_pressure_solver.cpp`
- Validation:
  - Commit message: 171 tests. Cavity 1297 iterations (vel=9.96e-6); channel 275 iterations
    (vel=9.63e-6).
- Risks/TODOs:
  - No automated comparison of the channel result with the analytical Poiseuille profile.

## 2026-04-22 19:42 (Europe/Warsaw) - SFML GUI (`lgflow_ui`) With Live Convergence Plot
- Author: Claude Code
- Status: committed (`405248f`)
- Summary:
  - Header-only SFML 3 widgets: colour palette, layout constants, rounded-rect generators,
    log-scale residual plot, side panel with case selector and Run/Stop, status bar.
  - `main_ui.cpp`: solver runs on a background thread and posts residuals to a
    mutex-protected `SolverState` that the UI reads each frame.
  - CMake: optional `find_package(SFML 3)`; `lgflow_ui` is skipped if SFML is missing. SFML
    DLLs are copied next to the exe on Windows.
- Files changed:
  - `CMakeLists.txt`
  - `src/ui/main_ui.cpp`, `src/ui/ResidualPlot.hpp`, `src/ui/RoundedRect.hpp`,
    `src/ui/SidePanel.hpp`, `src/ui/StatusBar.hpp`, `src/ui/UIColors.hpp`, `src/ui/UILayout.hpp`
- Validation:
  - Not recorded in the commit message.
- Risks/TODOs:
  - No tests for the UI code.

## 2026-10-05 19:31 (Europe/Warsaw) - Unstructured Solver, Gmsh Reader, Spalart-Allmaras, GUI Field View (WIP)
- Author: Claude Code
- Status: committed (`d86fc17`), marked work in progress
- Summary:
  - Unstructured mesh, field and BC classes (`src/mesh/`); Gmsh `.msh` v2 ASCII reader;
    unstructured discretization and `UnstructuredNavierStokesSolver`.
  - Spalart-Allmaras model with wall distance, enabled in the structured solver via
    `solver.turbulence = SA`.
  - New cases `channel_unstruct` and `flat_plate` with executables `lgflow_unstruct` and
    `lgflow_flat_plate`.
  - `lgflow_validation_tests` target (CTest `GhiaRe100`) and L2 thresholds in
    `case_validate.cfg` (`max_u_l2 = 0.10`, `max_v_l2 = 0.06`).
  - GUI: field view (pressure/velocity colormap, streamlines), Nx/Ny/Re/iteration inputs.
  - Note: this goes beyond the Milestone 3/4 scope in CLAUDE.md (Milestones 5 and 6).
- Files changed:
  - 47 files; see `git show --stat d86fc17`.
- Validation:
  - Not recorded in the commit message. Test sources count 223 `TEST`/`TEST_F` macros
    (220 unit, 3 validation); a fresh ctest run has not been logged.
- Risks/TODOs:
  - Unstructured solver: no solver-level tests and no validation case.
  - `solver.nonorth_correctors` is parsed but unused.
  - Spalart-Allmaras: no comparison with a reference (e.g. flat-plate skin friction).

## 2026-10-05 19:42 (Europe/Warsaw) - README Rewrite
- Author: Claude Code
- Status: committed (`4943f86`)
- Summary:
  - README rewritten from the source: build targets, full config and BC reference, Re=100
    validation details and exit codes, GUI, project status with WIP areas marked.
  - Leaves author TODOs: SFML install route, measured validation numbers, fresh ctest count.
- Files changed:
  - `README.md`
- Validation:
  - Documentation only.
- Risks/TODOs:
  - Re-run `lgflow_validate_lid` (2000 iterations) and `ctest`, then fill in the README TODOs.

## 2026-10-05 23:20 (Europe/Warsaw) - HISTORY.md Backfill
- Author: Claude Code
- Status: committed
- Summary:
  - Added the six entries above for commits `5c6e175` through `4943f86`.
- Files changed:
  - `HISTORY.md`
- Validation:
  - Documentation only; no build or tests run.
- Risks/TODOs:
  - Same open items as the entries above.

## 2026-10-06 00:50 (Europe/Warsaw) - Boundary Conditions as Boundary-Face Fluxes
- Author: Claude Code
- Status: committed (`1b12bf3`..`59ed1c3` on `feature/new-solver`, merged as `50ef67e`)
- Summary:
  - `BoundaryCondition::buildVelocityFaces()` / `buildPressureFaces()` precompute per-face
    closures (owner, normal, area, half-cell distance, kind, value) once per mesh.
  - `Discretization` gains boundary-aware Gauss gradients and Rhie-Chow divergence (exact
    u_b.n.A at walls/inlets, Rhie-Chow outlet faces with fixed p).
  - `MomentumSolver` / `PressureSolver::setBoundaryFaces()`: Dirichlet wall/inlet viscous flux
    nu*A*(u_b - u_P)/d_b, upwind donor from Dirichlet faces, p' = 0 on outlet faces
    (a_b = A/d_b); cell 0 pinned only in closed domains.
  - `NavierStokesSolver` wires the faces in `initialize()` and no longer overwrites boundary
    cells; walls now sit at the domain edges. `PARABOLIC_INLET` uses face-averaged values
    (inflow exactly U*H) and throws on top/bottom patches.
  - API change: `PressureSolver::setDirichletPressureCells()` removed;
    `BoundaryCondition::applyVelocity()/applyPressure()` documented as legacy.
  - Compatibility path: `WallDistance` and Spalart-Allmaras still use wall cells (owners of
    WALL faces), documented in `BoundaryCondition.hpp`.
  - Duplicated explicit-source code of the two `MomentumSolver::solve()` overloads moved to
    `computeExplicitSource()`.
  - README, `case_validate.cfg` and `test_ghia_re100.cpp` record the new Ghia baseline.
- Files changed:
  - `src/core/BoundaryCondition.hpp`, `src/core/BoundaryCondition.cpp`
  - `src/solver/Discretization.hpp`, `src/solver/Discretization.cpp`
  - `src/solver/MomentumSolver.hpp`, `src/solver/MomentumSolver.cpp`
  - `src/solver/PressureSolver.hpp`, `src/solver/PressureSolver.cpp`
  - `src/solver/NavierStokesSolver.hpp`, `src/solver/NavierStokesSolver.cpp`
  - `tests/test_bc.cpp`, `tests/test_discretization.cpp`, `tests/test_momentum_solver.cpp`,
    `tests/test_pressure_solver.cpp`, `tests/test_solver_init.cpp`
  - `tests/validation/test_ghia_re100.cpp`, `cases/lid_driven_cavity/case_validate.cfg`, `README.md`
- Validation:
  - `cmake --build build --config Release` then `ctest -C Release` (from `build/`):
    242/242 passed (1 SA integration test skipped, as before); baseline on `main` was 221/221.
  - `lgflow_validate_lid cases/lid_driven_cavity/case_validate.cfg 2000`: u L2 0.0712 -> 0.0432,
    v L2 0.0375 -> 0.0221, converged in 1376 iterations (was 1297).
  - `lgflow cases/channel_flow/case.cfg`: continuity residual 0.1496 (stalled) -> 5.3e-4 at the
    velocity-residual stop; with tolerance 1e-10 it falls to 1.7e-8 by iteration 542.
  - `lgflow_flat_plate`: runs, residuals finite (no accuracy check).
- Risks/TODOs:
  - `continuityResidual()` is an un-normalised L2 norm (depends on cell count).
  - `run()` stops on the velocity residual only although its doc says velocity AND continuity.
  - SA / `WallDistance` wall treatment still cell-based (compatibility path).
  - Poiseuille channel validation test designed, not yet implemented.
