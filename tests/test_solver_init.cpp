#include "solver/NavierStokesSolver.hpp"
#include "utils/Config.hpp"

#include <gtest/gtest.h>
#include <stdexcept>
#include <filesystem>
#include <fstream>

// ── Helpers ───────────────────────────────────────────────────────────────────

namespace
{
/// Returns a NavierStokesSolver constructed from an empty Config.
/// All solver parameters will use their built-in defaults, which are valid.
NavierStokesSolver makeDefaultSolver()
{
    Config cfg;  // empty: all keys will use constructor defaults
    return NavierStokesSolver(cfg);
}
}  // namespace

// ── Pre-initialization guards ─────────────────────────────────────────────────

TEST(NavierStokesSolverInit, Pressure_BeforeInit_Throws)
{
    auto solver = makeDefaultSolver();
    // pressure() must throw std::logic_error before initialize() is called
    EXPECT_THROW(static_cast<void>(solver.pressure()), std::logic_error);
}

TEST(NavierStokesSolverInit, Velocity_BeforeInit_Throws)
{
    auto solver = makeDefaultSolver();
    EXPECT_THROW(static_cast<void>(solver.velocity()), std::logic_error);
}

TEST(NavierStokesSolverInit, Step_BeforeInit_Throws)
{
    auto solver = makeDefaultSolver();
    EXPECT_THROW(solver.step(0.01), std::logic_error);
}

TEST(NavierStokesSolverInit, Run_BeforeInit_Throws)
{
    auto solver = makeDefaultSolver();
    EXPECT_THROW(solver.run(1), std::logic_error);
}

// ── Post-initialization access ────────────────────────────────────────────────

TEST(NavierStokesSolverInit, Initialize_DoesNotThrow)
{
    auto solver = makeDefaultSolver();
    EXPECT_NO_THROW(solver.initialize());
}

TEST(NavierStokesSolverInit, Pressure_AfterInit_DoesNotThrow)
{
    auto solver = makeDefaultSolver();
    solver.initialize();
    EXPECT_NO_THROW(static_cast<void>(solver.pressure()));
}

TEST(NavierStokesSolverInit, Velocity_AfterInit_DoesNotThrow)
{
    auto solver = makeDefaultSolver();
    solver.initialize();
    EXPECT_NO_THROW(static_cast<void>(solver.velocity()));
}

TEST(NavierStokesSolverInit, Step_AfterInit_DoesNotThrow)
{
    auto solver = makeDefaultSolver();
    solver.initialize();
    EXPECT_NO_THROW(solver.step(0.01));
}

TEST(NavierStokesSolverInit, Run_AfterInit_DoesNotThrow)
{
    auto solver = makeDefaultSolver();
    solver.initialize();
    // Run 1 stub iteration: should not throw
    EXPECT_NO_THROW(solver.run(1));
}

// ── Field sizing post-initialization ─────────────────────────────────────────

TEST(NavierStokesSolverInit, Pressure_AfterInit_SizeMatchesDefaultMesh)
{
    auto solver = makeDefaultSolver();
    solver.initialize();
    // Default mesh is 16x16 = 256 cells (from constructor defaults)
    EXPECT_EQ(solver.pressure().size(), 16 * 16);
}

TEST(NavierStokesSolverInit, Velocity_AfterInit_SizeMatchesDefaultMesh)
{
    auto solver = makeDefaultSolver();
    solver.initialize();
    EXPECT_EQ(solver.velocity().size(), 16 * 16);
}

// ── Re-initialization ─────────────────────────────────────────────────────────

TEST(NavierStokesSolverInit, ReInitialize_DoesNotThrow)
{
    // Calling initialize() twice should re-build the mesh and fields cleanly.
    auto solver = makeDefaultSolver();
    solver.initialize();
    EXPECT_NO_THROW(solver.initialize());
}

TEST(NavierStokesSolverInit, Residual_AlwaysAccessible)
{
    // residual() is valid at any lifecycle stage (returns current value, no throw)
    auto solver = makeDefaultSolver();
    EXPECT_NO_THROW(static_cast<void>(solver.residual()));
    solver.initialize();
    EXPECT_NO_THROW(static_cast<void>(solver.residual()));
}

// ── Invalid mesh config ───────────────────────────────────────────────────────

TEST(NavierStokesSolverInit, Initialize_ZeroNx_Throws)
{
    // Write a temporary config file with mesh.Nx = 0.
    // NavierStokesSolver::initialize() must propagate std::invalid_argument
    // from Mesh::load() when a non-positive dimension is supplied.
    namespace fs = std::filesystem;
    const fs::path tmpFile =
        fs::temp_directory_path() / "lgflow_test_zero_nx.cfg";

    {
        std::ofstream ofs(tmpFile);
        ASSERT_TRUE(ofs.is_open()) << "Cannot open temp config file for writing";
        ofs << "mesh.Nx = 0\n";
        ofs << "mesh.Ny = 16\n";
        ofs << "mesh.Lx = 1.0\n";
        ofs << "mesh.Ly = 1.0\n";
    }

    Config cfg;
    cfg.load(tmpFile.string());

    // Clean up before the assertion so we never leave a temp file behind.
    fs::remove(tmpFile);

    NavierStokesSolver solver(cfg);
    EXPECT_THROW(solver.initialize(), std::invalid_argument);
}

// ── BC wiring: config → solver → field ───────────────────────────────────────

/// After initialize(), the top boundary cells must carry the lid velocity
/// specified in the config (bc.top.value = 1.25 via fallback x-component).
TEST(NavierStokesSolverInit, BCWiring_TopLid_SetsTopRow)
{
    namespace fs = std::filesystem;
    const fs::path tmpFile =
        fs::temp_directory_path() / "lgflow_test_bc_wiring.cfg";

    {
        std::ofstream ofs(tmpFile);
        ASSERT_TRUE(ofs.is_open()) << "Cannot open temp config for writing";
        ofs << "bc.top.type    = WALL\n";
        ofs << "bc.top.value   = 1.25\n";   // fallback x-component
        ofs << "bc.bottom.type = WALL\n";
        ofs << "bc.left.type   = WALL\n";
        ofs << "bc.right.type  = WALL\n";
    }

    Config cfg;
    cfg.load(tmpFile.string());
    fs::remove(tmpFile);

    NavierStokesSolver solver(cfg);
    solver.initialize();

    const auto& vel = solver.velocity();
    const int Nx = 16, Ny = 16;  // defaults

    for (int i = 0; i < Nx; ++i) {
        EXPECT_NEAR(vel(i, Ny - 1).x(), 1.25, 1e-12) << "top row, i=" << i;
        EXPECT_NEAR(vel(i, Ny - 1).y(), 0.0,  1e-12) << "top row, i=" << i;
    }
}

/// An unknown BC type string must cause std::invalid_argument in the constructor.
TEST(NavierStokesSolverInit, BCWiring_InvalidType_Throws)
{
    namespace fs = std::filesystem;
    const fs::path tmpFile =
        fs::temp_directory_path() / "lgflow_test_bc_invalid.cfg";

    {
        std::ofstream ofs(tmpFile);
        ASSERT_TRUE(ofs.is_open()) << "Cannot open temp config for writing";
        ofs << "bc.top.type = NOTATYPE\n";
    }

    Config cfg;
    cfg.load(tmpFile.string());
    fs::remove(tmpFile);

    EXPECT_THROW(NavierStokesSolver solver(cfg), std::invalid_argument);
}
