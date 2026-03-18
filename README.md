# LG-Flow — 2D Incompressible Navier-Stokes CFD Solver

LG-Flow is an educational and research-grade CFD solver for 2D incompressible
viscous flow. It is built on the **Finite Volume Method (FVM)** with
**SIMPLE** (Semi-Implicit Method for Pressure-Linked Equations)
pressure-velocity coupling.

> **Status:** Skeleton — no physics implemented yet. All classes are stubs
> with documented interfaces. Build system and test infrastructure are in place.

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
- **Milestone 2** — Mesh & Fields: structured 2D mesh, Field arithmetic, VTK output
- **Milestone 3** — Laminar solver: SIMPLE loop, lid-driven cavity at Re=100
- **Milestone 4** — Validation: quantitative comparison against Ghia et al. (1982)
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
