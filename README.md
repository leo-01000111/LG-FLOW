# LG-Flow: 2D incompressible Navier-Stokes CFD solver

LG-Flow is a research and teaching code for 2D incompressible viscous flow. It
uses the finite volume method with SIMPLE (Semi-Implicit Method for
Pressure-Linked Equations) pressure-velocity coupling on a collocated grid:
velocity and pressure are both stored at cell centres, and Rhie-Chow face
interpolation is used in the pressure equation to avoid pressure-velocity
decoupling.

The structured-grid solver is the main code path and has unit tests and a
lid-driven cavity validation against Ghia et al. (1982) at Re = 100. The
unstructured-mesh solver and the Spalart-Allmaras turbulence model are
committed as work in progress (see [Project status](#project-status)).

## Dependencies

| Dependency | Version | Purpose |
|---|---|---|
| CMake | 3.20 or newer | Build system |
| C++ compiler | GCC 11 / Clang 14 / MSVC 19.29 or newer | C++20 support |
| [Eigen3](https://eigen.tuxfamily.org) | no minimum set in `find_package` | Linear algebra, sparse LU |
| [GoogleTest](https://github.com/google/googletest) | no minimum set in `find_package` | Unit tests (required by the build) |
| [SFML](https://www.sfml-dev.org) | 3.x | GUI (`lgflow_ui`), optional |

The compiler versions in the table are the ones the project was developed
against; CMake itself only requires C++20 support.

If SFML 3 is not found at configure time, the `lgflow_ui` target is skipped and
everything else still builds.

### Installing dependencies (Ubuntu/Debian)
```bash
sudo apt install cmake libeigen3-dev libgtest-dev
```

### Installing dependencies (macOS via Homebrew)
```bash
brew install cmake eigen googletest
```

### Installing dependencies (Windows, MSVC + vcpkg)
```powershell
# Core tools
winget install --id Kitware.CMake -e
winget install --id Microsoft.VisualStudio.2022.BuildTools -e

# Dependencies (Eigen3 + GoogleTest; add sfml:x64-windows for the GUI)
git clone https://github.com/microsoft/vcpkg C:\vcpkg
C:\vcpkg\bootstrap-vcpkg.bat
C:\vcpkg\vcpkg.exe install eigen3:x64-windows gtest:x64-windows
```

<!-- TODO: author to confirm the SFML 3 install route (vcpkg port name and version, Linux and macOS packages). -->

## Build instructions

```bash
# Configure (Debug build; adds ASan/UBSan on GCC and Clang)
cmake -B build -DCMAKE_BUILD_TYPE=Debug

# Configure (Release build; -O3 -march=native on GCC/Clang, /O2 /GL on MSVC)
cmake -B build-release -DCMAKE_BUILD_TYPE=Release

# Build
cmake --build build

# Run tests
cd build && ctest --output-on-failure

# Run the solver (from the repository root)
./build/lgflow cases/lid_driven_cavity/case.cfg
```

AddressSanitizer and UndefinedBehaviorSanitizer are only enabled for Debug
builds with GCC or Clang. MSVC builds use `/W4 /WX /permissive- /EHsc /utf-8`
instead and get no sanitizers. On GCC and Clang the warning flags are `-Wall
-Wextra -Wpedantic`.

### Build on Windows (PowerShell + vcpkg toolchain)
```powershell
cmake -S . -B build `
  -DCMAKE_BUILD_TYPE=Debug `
  -DCMAKE_TOOLCHAIN_FILE=C:/vcpkg/scripts/buildsystems/vcpkg.cmake `
  -DVCPKG_TARGET_TRIPLET=x64-windows

cmake --build build --config Debug
ctest --test-dir build --output-on-failure -C Debug
.\build\Debug\lgflow.exe cases\lid_driven_cavity\case.cfg
```

### Build targets

| Target | Source | Description |
|---|---|---|
| `flowcore_lib` | `src/core`, `src/solver`, `src/io`, `src/mesh`, `src/utils` | Static library shared by every executable below |
| `lgflow` | `src/main.cpp` | Structured-grid solver driven by a config file |
| `lgflow_validate_lid` | `src/validate/lid_driven_cavity_validate.cpp` | Lid-driven cavity run with comparison to Ghia et al. (1982) |
| `lgflow_unstruct` | `src/validate/unstructured_channel_validate.cpp` | Unstructured-mesh channel run (work in progress) |
| `lgflow_flat_plate` | `src/validate/flat_plate_sa_validate.cpp` | Flat plate run with the Spalart-Allmaras model (work in progress) |
| `lgflow_ui` | `src/ui/main_ui.cpp` | SFML 3 GUI, built only if SFML 3 is found |
| `lgflow_tests` | `tests/*.cpp` | GoogleTest unit tests |
| `lgflow_validation_tests` | `tests/validation/test_ghia_re100.cpp` | Ghia Re = 100 regression test, registered in CTest as `GhiaRe100` |
| `lgflow_poiseuille_validation` | `tests/validation/test_poiseuille_channel.cpp` | Poiseuille channel grid-refinement study, registered in CTest as `PoiseuilleChannel` (slow) |

All commands that read `cases/...` should be run from the repository root. CTest
sets the working directory for `GhiaRe100` and `PoiseuilleChannel` itself.

## Running the solver

```bash
lgflow <config-file>
```

`lgflow` loads the config, runs up to `solver.max_iter` SIMPLE iterations
(1000 if the key is absent) and stops early when the convergence criterion below
is met. It returns 1 on a missing config file or any exception, and 0
otherwise. Reaching `solver.max_iter` without converging only logs a warning and
still returns 0.

Files written to `output.dir` (default `output`):

| File | Description |
|---|---|
| `history.csv` | One row per iteration with columns `iter,vel_residual,cont_residual,pressure_residual` |
| `iter_<N>.vtu` | VTK snapshot every `output.vtk_interval` iterations |
| `result.vtu` | Final fields, written by `lgflow` after the run |

The `.vtu` files are ASCII XML (UnstructuredGrid), one vertex per cell centre,
with cell data `pressure` and `velocity`. They open in ParaView.

### Convergence criterion

The run stops only when both residuals are below their tolerances
(`NavierStokesSolver::isConverged()`, also used by the GUI):

- Velocity residual `||u_new - u_old|| / max(||u_old||, 1e-12)` below
  `solver.tolerance`.
- Relative continuity residual `R_k / R_0` below `solver.tolerance_continuity`.
  `R_k` is the L2 norm of the Rhie-Chow divergence of the predicted velocity
  (before the pressure correction), with exact fluxes on wall and inlet faces.
  `R_0` is the largest `R_k` of the first 5 iterations; until then the value is
  reported as 1. Being a ratio of norms over the same cells, it does not grow
  with the cell count, so it is comparable across mesh refinements.

Both values are written to `history.csv` and the log.

## Algorithm summary (structured solver)

- Grid: uniform Cartesian, `mesh.Nx` by `mesh.Ny` cells over `mesh.Lx` by
  `mesh.Ly`. Cell index is `i * Ny + j`. Velocity and pressure are both at cell
  centres.
- Momentum predictor: implicit viscous term, explicit convection and pressure
  gradient, sparse LU solve (Eigen). Velocity under-relaxation `alpha_u` is
  built into the matrix diagonal and right-hand side.
- Convection: `upwind` (first-order donor cell) or `central` (face average).
- Pressure correction: Poisson equation for `p'` with Rhie-Chow face
  velocities in the divergence, solved with sparse LU.
  `solver.pressure_corrections_per_step` repeats this solve within one SIMPLE
  step. The pressure is updated with `alpha_p`.
- Boundary conditions are applied as boundary-face fluxes: walls and inlets
  prescribe the velocity on the boundary face (half a cell from the first cell
  centre), outlets fix the pressure on the outlet face. No cell values are
  overwritten, so the walls sit exactly at the domain edges.
- Pressure reference: p' = 0 on OUTLET faces. If the case has no OUTLET,
  cell 0 is pinned.
- The time step is clamped each iteration to
  `min(dt, max_cfl_conv * h / max|u|, max_cfl_diff * h^2 / nu)`, with `h` the
  smaller cell size. The effective value is shown as `dt_eff` in the log.

## Configuration reference

Config files are plain text, one `key = value` per line. Lines starting with `#`
are comments. Keys not listed below are ignored.

### Mesh and solver (structured grid, `lgflow`)

| Key | Default | Description |
|---|---|---|
| `mesh.Nx`, `mesh.Ny` | `16`, `16` | Cells in x and y |
| `mesh.Lx`, `mesh.Ly` | `1.0`, `1.0` | Domain size |
| `solver.dt` | `0.01` | Time step (must be > 0) |
| `solver.max_iter` | `1000` in `lgflow` | Maximum SIMPLE iterations |
| `solver.tolerance` | `1e-6` | Velocity residual convergence threshold (must be > 0) |
| `solver.tolerance_continuity` | `1e-4` | Relative continuity residual threshold (must be > 0); both thresholds must be met to stop |
| `solver.rho` | `1.0` | Density (must be > 0) |
| `solver.nu` | `0.01` | Kinematic viscosity (must be >= 0) |
| `solver.alpha_u` | `0.7` | Velocity under-relaxation, in (0, 1] |
| `solver.alpha_p` | `0.3` | Pressure under-relaxation, in (0, 1] |
| `solver.convection_scheme` | `upwind` | `upwind` or `central` (case-insensitive) |
| `solver.max_cfl_conv` | `0.5` | Convective CFL limit used to clamp dt (must be > 0) |
| `solver.max_cfl_diff` | `0.5` | Diffusive CFL limit used to clamp dt (must be > 0) |
| `solver.pressure_corrections_per_step` | `1` | Pressure Poisson solves per SIMPLE step (must be >= 1) |
| `solver.turbulence` | `LAMINAR` | `SA` (or `SPALART-ALLMARAS`) enables the Spalart-Allmaras model; any other value runs laminar |
| `solver.sa.alpha_nu` | `0.7` | Under-relaxation for the SA variable, in (0, 1]; read only when SA is on |
| `solver.sa.nu_tilde_freestream` | `3 * solver.nu` | Initial and freestream value of the SA variable; read only when SA is on |

### Output

| Key | Default | Description |
|---|---|---|
| `output.dir` | `output` | Directory for `history.csv` and `.vtu` files |
| `output.vtk_interval` | `100` | Write `iter_<N>.vtu` every N iterations (must be > 0) |

### Boundary conditions

Boundaries are named `left`, `right`, `bottom` and `top`. Each is configured
with:

```
bc.<left|right|bottom|top>.type    = WALL | INLET | OUTLET | SYMMETRY | PARABOLIC_INLET
bc.<side>.value                    = <x-component>     # shorthand for value_x
bc.<side>.value_x                  = <x-component>
bc.<side>.value_y                  = <y-component>
```

A side without a `.type` key gets no boundary treatment. Type names are
case-insensitive. `value_x` takes precedence over `value`; both default to 0.

| Type | Velocity | Pressure |
|---|---|---|
| `WALL` | Dirichlet on the boundary face: `(value_x, value_y)`, so a non-zero value gives a moving wall (the lid in the cavity) | zero gradient |
| `INLET` | Dirichlet on the boundary face: `(value_x, value_y)` | zero gradient |
| `PARABOLIC_INLET` (`left` or `right` only) | Poiseuille profile `6 U y (Ly - y) / Ly^2` with `U = value_x`, averaged over each face so the inflow is exactly `U * Ly` | zero gradient |
| `OUTLET` | zero gradient | fixed to 0 on the outlet face |
| `SYMMETRY` | zero gradient | zero gradient |

`SYMMETRY` leaves the boundary values unchanged. It does not zero the normal
velocity component, so it acts as a zero-gradient (slip-like) boundary only.

Each boundary face belongs to exactly one side, so there are no corner
conflicts. The Spalart-Allmaras model still treats the cells next to `WALL`
faces as its wall cells (wall distance to their centres, `nu_tilde = 0` there);
this compatibility path is documented in `BoundaryCondition.hpp`.

### Unstructured solver (`lgflow_unstruct`)

| Key | Default | Description |
|---|---|---|
| `mesh.file` | `mesh.msh` | Path to a Gmsh `.msh` file (format 2, ASCII) |
| `solver.dt`, `solver.rho`, `solver.nu`, `solver.tolerance`, `solver.alpha_u`, `solver.alpha_p` | as above | Same meaning as in the structured solver |
| `solver.nonorth_correctors` | `2` | Read into a member variable; no other use of it was found in the source |
| `output.dir`, `output.vtk_interval` | as above | Same meaning as above |
| `bc.<patch>.type`, `.value`, `.value_x`, `.value_y` | none | `<patch>` is the Gmsh physical group name (for example `inlet`, `wall_top`). Types: `WALL`, `INLET`, `OUTLET`, `SYMMETRY` |

On the unstructured mesh only `WALL` and `INLET` set velocity (Dirichlet);
`OUTLET` and `SYMMETRY` leave it unchanged, and `OUTLET` cells are the pressure
reference cells (cell 0 if there is no OUTLET). `PARABOLIC_INLET` is accepted by
the config parser but has no effect there.

The unstructured solver always uses explicit upwind convection and does not read
`solver.convection_scheme`, the CFL keys, `solver.pressure_corrections_per_step`
or the turbulence keys.

## Validation: lid-driven cavity at Re = 100 and Re = 1000

The validation benchmark is the lid-driven cavity at Re = 100 (`nu = 0.01`,
lid speed 1, unit square), compared with the centreline data in
`cases/lid_driven_cavity/ghia1982_re100.csv`:

> Ghia, U., Ghia, K.N., Shin, C.T. (1982). *High-Re solutions for
> incompressible flow using the Navier-Stokes equations and a multigrid
> method.* Journal of Computational Physics, 48(3), 387-411.

The Re = 1000 case (`case_validate_re1000.cfg`, `nu = 0.001`) uses
`cases/lid_driven_cavity/ghia1982_re1000.csv` from the same paper (Tables I and
II). The config key `validation.reference` selects the reference file. The
default is the Re = 100 file.

### Running the validation executable

Run from the repository root so the relative paths to `cases/` resolve:

```powershell
# Default: cases/lid_driven_cavity/case_validate.cfg, at most 500 iterations
.\build\Debug\lgflow_validate_lid.exe

# Custom config and/or iteration count
.\build\Debug\lgflow_validate_lid.exe cases\lid_driven_cavity\case_validate.cfg 2000
```

The default iteration cap is `min(solver.max_iter, 500)`. `case_validate.cfg`
sets `solver.max_iter = 2000`, which is also the count the `GhiaRe100` CTest
uses, so pass `2000` explicitly to match it.

`case_validate.cfg` uses a 32 x 32 mesh, upwind convection, `dt = 0.01`,
`alpha_u = 0.7` and `alpha_p = 0.3`.

Files written to `output/lid_driven_cavity/`:

| File | Description |
|---|---|
| `history.csv` | Per-iteration residuals (`iter,vel_residual,cont_residual,pressure_residual`) |
| `centerline.csv` | Sampled u along the vertical centreline and v along the horizontal centreline, columns `axis,coord,value` |
| `validation_metrics.csv` | `metric,value` rows: `ref_available`, `final_vel_residual`, `final_cont_residual`, `u_l2`, `u_linf`, `v_l2`, `v_linf`, `iterations` |
| `iter_*.vtu` | VTK snapshots every `output.vtk_interval` iterations |

Each reference point is matched to the nearest sampled coordinate. The
executable computes both the L2 error (root mean square) and the L-infinity
error (maximum absolute) for u and v, prints them, and writes them to
`validation_metrics.csv`. Only the L2 errors and the continuity residual can be
checked against thresholds. L-infinity is reported and never compared with a
limit. If the reference CSV is missing, the four error rows are written as
`nan`.

### Pass/fail thresholds

| Key | Description |
|---|---|
| `validation.reference` | Reference CSV to compare with (default `cases/lid_driven_cavity/ghia1982_re100.csv`). |
| `validation.max_u_l2` | Fails if the u-centreline L2 error is larger. Ignored when the reference CSV is missing. |
| `validation.max_v_l2` | Fails if the v-centreline L2 error is larger. Ignored when the reference CSV is missing. |
| `validation.max_cont_residual` | Fails if the final relative continuity residual is larger. |

`case_validate.cfg` sets `validation.max_u_l2 = 0.10` and
`validation.max_v_l2 = 0.06`. Comments in that file and in
`tests/validation/test_ghia_re100.cpp` record measured errors of 0.043 (u) and
0.022 (v) for the 32 x 32 upwind setup (face-based boundary conditions;
the earlier cell-overwrite boundaries gave 0.071 and 0.037).

Measured on 2026-10-27 (Release, MSVC) with
`lgflow_validate_lid cases/lid_driven_cavity/case_validate.cfg 2000`:

| Metric | Value |
|---|---|
| Iterations to convergence | 1376 |
| Final velocity residual | 9.99e-6 |
| Final relative continuity residual | 7.15e-6 |
| u centreline L2 / Linf | 0.0432 / 0.113 |
| v centreline L2 / Linf | 0.0221 / 0.0389 |

Exit codes of `lgflow_validate_lid`:

| Code | Meaning |
|---|---|
| `0` | Run completed, and no threshold was set or none was exceeded |
| `1` | Fatal error (config missing, solver exception, file write failure) |
| `2` | At least one threshold was exceeded (`FAIL` is printed with the values) |

When thresholds are present and all pass, the executable still prints a line
saying the metrics are informational. The exit code is what indicates the
result. The other executables (`lgflow`, `lgflow_unstruct`, `lgflow_flat_plate`)
only return 0 or 1 and never fail on accuracy.

`lgflow_validation_tests` (CTest name `GhiaRe100`) runs the same case and has
three tests: the velocity residual must fall below `1e-4`, and the u and v
centreline L2 errors must be below 0.10 and 0.06. It is skipped if the case
files are not found.

### Re = 1000

CTest `GhiaRe1000` (labels `slow` and `validation`, about 3 minutes) runs
`lgflow_validate_lid cases/lid_driven_cavity/case_validate_re1000.cfg 20000`
and fails on a non-zero exit code. The config uses a 64 x 64 mesh with upwind
convection and sets `validation.max_u_l2 = 0.15`, `validation.max_v_l2 = 0.18`
and `validation.max_cont_residual = 1e-4`.

Measured on 2026-10-27 (Release, MSVC), upwind convection:

| Mesh | u L2 / Linf | v L2 / Linf | Iterations | Time |
|---|---|---|---|---|
| 32 x 32 | 0.172 / 0.334 | 0.212 / 0.371 | 4402 | 35 s |
| 64 x 64 | 0.103 / 0.184 | 0.124 / 0.201 | 5111 | 200 s |
| 128 x 128 | 0.063 / 0.110 | 0.070 / 0.100 | 9616 | 2371 s |

The errors are several times larger than at Re = 100 and fall by only about
1.7 per mesh doubling. The numerical diffusion of first-order upwind
convection is large compared with the physical viscosity at Re = 1000. Central
differencing diverges on the 64 x 64 mesh, because the cell Péclet number
U·dx/ν = 16 is far above the limit of 2 for a bounded central scheme. Closer
agreement needs a bounded higher-order convection scheme, which is not
implemented.

## Validation: Poiseuille channel flow (grid refinement)

`lgflow_poiseuille_validation` (CTest name `PoiseuilleChannel`, labels `slow`
and `validation`) solves `cases/channel_flow/case_validate.cfg`: a uniform
inflow U = 1 m/s between two no-slip walls at Re_H = U·H/ν = 100. It checks that
the flow develops into the analytical plane Poiseuille profile
u(y) = 6·U·(y/H)·(1 − y/H) (White 2006, §6-2).

The domain length is 3 times the parallel-plate entrance length of Chen (1973),
rounded up to whole coarse cells (Lx = 14 H). The test runs the case on meshes
with Ny = 8, 16, 32, … cells across the channel and the same cell aspect ratio
dx/dy = 4 on every level. It checks the following on each level:

- the solver converges (velocity residual < 1e-8, relative continuity < 1e-6);
- the flow rate Q = Σ u·dy, summed from cell-centre velocities, equals U·H
  within `validation.max_q_error` at mid-channel, the profile station and the
  last cell column;
- the Q error at x = 1 H and the RMSE of u(y) against the analytical profile at
  x = 0.9 Lx (divided by u_max = 1.5 U) are both lower than on the previous
  level.

Cell-centre velocities match the conserved face fluxes only up to a term of
order dx²·∂²u/∂x². That term is large near the inlet, where the uniform inflow
meets the walls, so x = 1 H only has to converge under refinement. It has to
meet the tolerance on the finest level only.

After the last level, the observed order of accuracy of the two finest levels,
p = log2(E_coarse / E_fine), must lie between `validation.min_order` and
`validation.max_order`. The finest RMSE must be below `validation.max_rel_rmse`,
and the finest Q error at x = 1 H below `validation.max_q_error`. The first three levels always run. A further level
runs only if its predicted time fits the budget set by the environment variable
`LGFLOW_VALIDATION_BUDGET_S` (default 1800 s). The CTest timeout is 2400 s, so a
larger budget also needs `ctest --timeout <s>`.

```bash
ctest --test-dir build -C Release -R PoiseuilleChannel --output-on-failure
```

Measured on 2026-10-27 (Release, MSVC, default budget). Q error columns are
|Q − U·H| / (U·H):

| Ny | Nx | Time [s] | Relative RMSE | Q error at 1 H | Q error, developed stations (max) |
|---|---|---|---|---|---|
| 8 | 28 | 0.5 | 1.02e-2 | 6.4e-3 | 2.1e-6 |
| 16 | 56 | 4.4 | 2.65e-3 | 5.1e-4 | 2.2e-7 |
| 32 | 112 | 38 | 6.79e-4 | 1.2e-5 | 2.6e-8 |
| 64 | 224 | 460 | 1.80e-4 | 1.2e-6 | 3.1e-9 |

The observed order is 1.94, 1.96 and 1.92 for the successive pairs. It is
second order, as expected from the central-differenced diffusion term, which
dominates in fully developed flow, where convection vanishes. The default budget
runs four levels, about 8.5 minutes in total. A fifth level (Ny = 128) would take
an estimated 1.5 hours. To skip the test in a quick run, use
`ctest -LE slow`.

## Test suite

`lgflow_tests` covers the mesh, field algebra, discretisation operators,
boundary conditions, config parsing, logger, momentum and pressure solvers,
solver parameter checks and run loop, convection schemes, VTK writer,
unstructured mesh and field containers, wall distance and the Spalart-Allmaras
closure functions.

A Release `ctest` run on 2026-10-27, excluding the slow `PoiseuilleChannel`
(`ctest -C Release -E PoiseuilleChannel`), passed 250 of 250 tests in 198 s,
171 s of it `GhiaRe1000`. One test,
`SpallartAllmaras.NavierStokesSolver_Integration_OneStep`, is skipped.
`ctest -LE slow` skips both slow validation cases.

There are no tests that run the unstructured solver or the Gmsh reader on a full
case.

## GUI (`lgflow_ui`)

`lgflow_ui` is an SFML 3 front end for the structured solver. It is built only if
`find_package(SFML 3 ...)` succeeds. The solver runs on a worker thread and the
window shows a live convergence plot (velocity and continuity residuals) and the
field, refreshed every 25 iterations, with a pressure or velocity colormap and a
streamline display mode. The side panel selects one of two cases, Lid-Driven
Cavity Re = 100 (`cases/lid_driven_cavity/case_validate.cfg`) and Channel Flow
Re = 100 (`cases/channel_flow/case.cfg`), and lets you override `Nx`, `Ny`, Re
and the iteration limit. It looks for `cases/` in the working directory and then
in the source directory compiled into the binary. It does not run the
unstructured solver or the turbulence cases.

## Example cases

| Case | Solver | Description |
|---|---|---|
| `cases/lid_driven_cavity/case.cfg` | structured | 64 x 64 cavity, Re = 100, `max_iter = 5000`, `tolerance = 1e-6` |
| `cases/lid_driven_cavity/case_validate.cfg` | structured | 32 x 32 cavity with validation thresholds |
| `cases/channel_flow/case.cfg` | structured | 64 x 16 plane channel, Re = 100 (based on mean velocity and height), `PARABOLIC_INLET` on the left, `OUTLET` on the right |
| `cases/flat_plate/case.cfg` | structured, `solver.turbulence = SA` | 64 x 32 cells over 1.0 by 0.2, `INLET` left, `OUTLET` right, `WALL` bottom, `SYMMETRY` top |
| `cases/channel_unstruct/case.cfg` | unstructured | Channel on the Gmsh mesh `cases/channel_unstruct/channel.msh`, patches `inlet`, `outlet`, `wall_bottom`, `wall_top` |

`case.cfg` starts from the analytical Poiseuille profile at the inlet. The
automated comparison uses `case_validate.cfg` instead, which has a uniform inlet
(see "Validation: Poiseuille channel flow").

## Project structure

Regenerated from `git ls-files` (the repository is `LG-FLOW`):

```
LG-FLOW/
├── CMakeLists.txt
├── HISTORY.md                 # engineering log (last entry covers Phase 9)
├── README.md
├── cases/
│   ├── channel_flow/          # case.cfg, case_validate.cfg
│   ├── channel_unstruct/      # case.cfg, channel.msh
│   ├── flat_plate/            # case.cfg
│   └── lid_driven_cavity/     # case.cfg, case_validate.cfg, ghia1982_re100.csv
├── src/
│   ├── main.cpp               # lgflow
│   ├── core/                  # Mesh, Field, BoundaryCondition (structured)
│   ├── mesh/                  # UnstructuredMesh, UnstructuredField, UnstructuredBoundaryCondition
│   ├── solver/                # NavierStokesSolver, MomentumSolver, PressureSolver, Discretization,
│   │                          #   UnstructuredNavierStokesSolver, UnstructuredDiscretization,
│   │                          #   SpallartAllmaras, WallDistance
│   ├── io/                    # MeshReader (stub), GmshReader, VTKWriter
│   ├── utils/                 # Logger, Config
│   ├── validate/              # sources of lgflow_validate_lid, lgflow_unstruct, lgflow_flat_plate
│   └── ui/                    # SFML 3 GUI (main_ui.cpp and header-only widgets)
└── tests/
    ├── CMakeLists.txt
    ├── test_*.cpp             # GoogleTest unit tests (lgflow_tests)
    └── validation/            # test_ghia_re100.cpp (lgflow_validation_tests)
```

`io/MeshReader` is still a stub that returns `false`; structured meshes are
generated from the `mesh.*` keys instead.

## Project status

| Area | Status |
|---|---|
| Structured mesh, fields, discretisation operators, BCs | Implemented and unit tested |
| SIMPLE loop with implicit viscous term, Rhie-Chow, CFL clamp, upwind and central convection | Implemented and unit tested |
| Lid-driven cavity vs Ghia et al. (1982), Re = 100 and Re = 1000 | Implemented, with thresholds on the L2 errors enforced by `lgflow_validate_lid`, `GhiaRe100` and `GhiaRe1000` (Re = 1000 accuracy limited by upwind convection) |
| INLET, OUTLET, PARABOLIC_INLET with pressure reference cells; `channel_flow` case | Implemented. Developed profile vs analytical Poiseuille, second-order grid convergence enforced by `PoiseuilleChannel` |
| GUI (`lgflow_ui`) | Implemented for the two structured cases listed above |
| Unstructured-mesh solver and Gmsh reader | Work in progress. Mesh and field containers have unit tests. The solver has no tests and no validation case. |
| Spalart-Allmaras model with wall distance | Work in progress. Implemented in the structured solver only. The closure functions, wall distance and a one-step smoke test on an 8 x 8 mesh are tested. No result has been compared with a reference such as a flat plate skin-friction correlation. |
| SYMMETRY boundary | Zero-gradient only (see the boundary table) |

<!-- TODO: author to add measured results for channel_unstruct and flat_plate once they have been run and checked against a reference. -->

`HISTORY.md` is the engineering log, with one entry per change.

## License

Proprietary / Research Use Only. Contact the author for licensing details.
