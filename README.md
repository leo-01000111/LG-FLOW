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

All commands that read `cases/...` should be run from the repository root. CTest
sets the working directory for `GhiaRe100` itself.

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

The run stops when the velocity residual drops below `solver.tolerance`. The
velocity residual is `||u_new - u_old|| / max(||u_old||, 1e-12)`, measured on
the committed velocity after boundary conditions are applied. The continuity
residual is written to `history.csv` and the log as a diagnostic and does not
stop the run. It is the norm of the Rhie-Chow divergence of the predicted
velocity, taken before the pressure correction. For cases with INLET/OUTLET
boundaries, a comment in `NavierStokesSolver.cpp` notes that it keeps a non-zero
offset on the collocated grid even at steady state, which is why it is not used
as a stopping test.

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
- Boundary conditions are applied to the predicted velocity before the pressure
  correction and again after it.
- Pressure reference: p' is pinned to zero on all OUTLET boundary cells. If the
  case has no OUTLET, cell 0 is used.
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
| `WALL` | Dirichlet: the boundary cells are set to `(value_x, value_y)`, so a non-zero value gives a moving wall (the lid in the cavity) | zero gradient |
| `INLET` | Dirichlet: boundary cells set to `(value_x, value_y)` | zero gradient |
| `PARABOLIC_INLET` (`left` or `right` only) | Poiseuille profile `6 U y (Ly - y) / Ly^2` with `U = value_x` | zero gradient |
| `OUTLET` | zero gradient | boundary cells fixed to 0, and used as the pressure-correction reference cells |
| `SYMMETRY` | zero gradient | zero gradient |

`SYMMETRY` leaves the boundary values unchanged. It does not zero the normal
velocity component, so it acts as a zero-gradient (slip-like) boundary only.

Corner cells are written in the order left, right, bottom, top, so the top
boundary wins at the top corners.

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

## Validation: lid-driven cavity at Re = 100

The validation benchmark is the lid-driven cavity at Re = 100 (`nu = 0.01`,
lid speed 1, unit square), compared with the centreline data in
`cases/lid_driven_cavity/ghia1982_re100.csv`:

> Ghia, U., Ghia, K.N., Shin, C.T. (1982). *High-Re solutions for
> incompressible flow using the Navier-Stokes equations and a multigrid
> method.* Journal of Computational Physics, 48(3), 387-411.

Only the Re = 100 reference data is in the repository. No Re = 1000 comparison
is implemented.

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
| `validation.max_u_l2` | Fails if the u-centreline L2 error is larger. Ignored when the reference CSV is missing. |
| `validation.max_v_l2` | Fails if the v-centreline L2 error is larger. Ignored when the reference CSV is missing. |
| `validation.max_cont_residual` | Fails if the final continuity residual is larger. |

`case_validate.cfg` sets `validation.max_u_l2 = 0.10` and
`validation.max_v_l2 = 0.06`. Comments in that file and in
`tests/validation/test_ghia_re100.cpp` record measured errors of 0.071 (u) and
0.037 (v) for the 32 x 32 upwind setup.

<!-- TODO: author to re-run lgflow_validate_lid on the committed code and enter the measured u_l2, u_linf, v_l2, v_linf and iteration count here. -->

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

## Test suite

`lgflow_tests` covers the mesh, field algebra, discretisation operators,
boundary conditions, config parsing, logger, momentum and pressure solvers,
solver parameter checks and run loop, convection schemes, VTK writer,
unstructured mesh and field containers, wall distance and the Spalart-Allmaras
closure functions. Counting `TEST` and `TEST_F` macros gives 223 test cases
across `tests/` (220 in `lgflow_tests` and 3 in `lgflow_validation_tests`).

There are no tests that run the unstructured solver or the Gmsh reader on a full
case.

<!-- TODO: author to confirm the pass count from a fresh ctest run (the 223 above is a source count, not a run result). -->

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

The channel_flow case file names the analytical Poiseuille profile as its
reference, but no executable compares the result with it.

## Project structure

Regenerated from `git ls-files` (the repository is `LG-FLOW`):

```
LG-FLOW/
├── CMakeLists.txt
├── HISTORY.md                 # engineering log (last entry covers Phase 9)
├── README.md
├── cases/
│   ├── channel_flow/          # case.cfg
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
| Lid-driven cavity vs Ghia et al. (1982), Re = 100 | Implemented, with thresholds on the L2 errors enforced by `lgflow_validate_lid` and `GhiaRe100` |
| INLET, OUTLET, PARABOLIC_INLET with pressure reference cells; `channel_flow` case | Implemented; no comparison with the analytical profile is automated |
| GUI (`lgflow_ui`) | Implemented for the two structured cases listed above |
| Unstructured-mesh solver and Gmsh reader | Work in progress. Mesh and field containers have unit tests. The solver has no tests and no validation case. |
| Spalart-Allmaras model with wall distance | Work in progress. Implemented in the structured solver only. The closure functions, wall distance and a one-step smoke test on an 8 x 8 mesh are tested. No result has been compared with a reference such as a flat plate skin-friction correlation. |
| SYMMETRY boundary | Zero-gradient only (see the boundary table) |

<!-- TODO: author to add measured results for channel_unstruct and flat_plate once they have been run and checked against a reference. -->

`HISTORY.md` records work up to Phase 9. The later commits (SIMPLE fixes,
INLET/OUTLET, GUI, unstructured solver, Spalart-Allmaras) have no entries there.

## License

Proprietary / Research Use Only. Contact the author for licensing details.
