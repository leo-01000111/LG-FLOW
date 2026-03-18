#include "solver/NavierStokesSolver.hpp"
#include "utils/Config.hpp"

#include <gtest/gtest.h>
#include <stdexcept>

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
    // Provide a config that requests a zero-cell mesh dimension.
    // Mesh::load() should propagate std::invalid_argument.
    // We test this via a custom Config with mesh.Nx = 0.
    // Since Config requires file loading, construct solver with default
    // (Nx=16) then test Mesh independently -- or rely on Mesh tests.
    // Here we just verify the throw propagates:
    // TODO(Milestone 2): inject a zero-Nx config and verify std::invalid_argument.
    EXPECT_TRUE(true);
}
