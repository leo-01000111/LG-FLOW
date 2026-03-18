# LG-Flow — 2D Incompressible Navier-Stokes CFD Solver

LG-Flow is an educational and research-grade CFD solver for 2D incompressible
viscous flow. It is built on the **Finite Volume Method (FVM)** with
**SIMPLE** (Semi-Implicit Method for Pressure-Linked Equations)
pressure-velocity coupling.

> **Status:** Phase 8 — upwind convection scheme, CFL safety clamp, and
> improved validation harness (validation_metrics.csv).

---

## Dependencies

| Dependency | Minimum version | Purpose |
|---|---|---|
| CMake | 3.20 | Build system |
| C++ compiler | GCC 11 / Clang 14 / MSVC 19.29 | C++20 support |
| [Eigen3](https://eigen.tuxfamily.org) | 3.3 | Linear algebra |
| [GoogleTest](https://github.com/google/googletest) | 1.12 | Unit testing |

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

# Dependencies (Eigen3 + GoogleTest)
git clone https://github.com/microsoft/vcpkg C:\vcpkg
C:\vcpkg\bootstrap-vcpkg.bat
C:\vcpkg\vcpkg.exe install eigen3:x64-windows gtest:x64-windows
```

---

## Build instructions

```bash
# Configure (Debug build with ASan/UBSan)
cmake -B build -DCMAKE_BUILD_TYPE=Debug

# Configure (Release build with -O3)
cmake -B build-release -DCMAKE_BUILD_TYPE=Release

# Build
cmake --build build

# Run tests
cd build && ctest --output-on-failure

# Run solver
./build/lgflow cases/lid_driven_cavity/case.cfg
```

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

### Running the validation executable (Windows)

Run from the `FlowCore/` directory so relative paths to `cases/` resolve correctly:

```powershell
# Default: loads cases/lid_driven_cavity/case_validate.cfg, caps at 500 iterations
.\build\Debug\lgflow_validate_lid.exe

# Custom config and/or iteration count
.\build\Debug\lgflow_validate_lid.exe cases\lid_driven_cavity\case_validate.cfg 1000
```

**Artifacts produced** (relative to `FlowCore/`):

| Artifact | Description |
|---|---|
| `output/lid_driven_cavity/history.csv` | Per-iteration residuals (iter, vel, cont, pressure) |
| `output/lid_driven_cavity/centerline.csv` | Sampled u(x=Lx/2,y) and v(x,y=Ly/2) with columns `axis,coord,value` |
| `output/lid_driven_cavity/validation_metrics.csv` | Machine-readable summary: final residuals, L2/Linf errors, iteration count |
| `output/lid_driven_cavity/iter_*.vtu` | VTK snapshots every `output.vtk_interval` iterations, openable in ParaView |

**Benchmark comparison mode:**

If `cases/lid_driven_cavity/ghia1982_re100.csv` is present (included in repo),
the executable automatically loads it and reports L2 and Linf error metrics
against Ghia et al. (1982) Re=100 data. Metrics are informational only.

**New config keys (Phase 8):**

| Key | Default | Description |
|---|---|---|
| `solver.convection_scheme` | `upwind` | Convection discretization: `upwind` (donor-cell, 1st order, stable) or `central` (Gauss, 2nd order) |
| `solver.max_cfl_conv` | `0.5` | Max convective CFL number; clamps dt to `max_cfl_conv * h / ||u||_max` |
| `solver.max_cfl_diff` | `0.5` | Max diffusive CFL number; clamps dt to `max_cfl_diff * h² / ν` |

---

## Project structure

```
FlowCore/
├── src/
│   ├── core/          # Mesh, Field, BoundaryCondition
│   ├── solver/        # NavierStokesSolver, PressureSolver, Discretization
│   ├── io/            # MeshReader, VTKWriter
│   └── utils/         # Logger, Config
├── tests/             # GoogleTest unit tests
└── cases/             # Example simulation configurations
    └── lid_driven_cavity/
```

---

## Planned features

- **Milestone 1** ✅ — Skeleton: project structure, CMake, stub classes
- **Milestone 2** ✅ — Mesh & Fields: structured 2D mesh, Field arithmetic, VTK output
- **Milestone 3** ✅ — Laminar solver: SIMPLE loop, residual tracking
- **Milestone 4** 🔄 — Validation harness: upwind convection, CFL clamp, centerline CSV + Ghia et al. (1982) comparison
- **Milestone 5** — Unstructured mesh support *(future)*
- **Milestone 6** — Turbulence modeling: Spalart-Allmaras *(future)*

---

## Canonical validation case

The **lid-driven cavity** at Re=100 and Re=1000 is used as the primary
validation benchmark, compared against the reference data from:

> Ghia, U., Ghia, K.N., Shin, C.T. (1982). *High-Re solutions for
> incompressible flow using the Navier-Stokes equations and a multigrid
> method.* Journal of Computational Physics, 48(3), 387–411.

---

## License

Proprietary / Research Use Only — contact author for licensing details.
