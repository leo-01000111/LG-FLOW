#include "solver/NavierStokesSolver.hpp"
#include "utils/Config.hpp"

#include <filesystem>
#include <fstream>
#include <gtest/gtest.h>
#include <stdexcept>
#include <string>

// Tests that NavierStokesSolver validates solver/output parameters at
// construction time, and validates dt in step() and maxIter in run().

// ── Helper ────────────────────────────────────────────────────────────────────

namespace
{

/// Writes a single-line config file to a temp path, loads it, removes the
/// file, and returns the Config.  Caller receives only the loaded Config.
Config singleKeyConfig(const std::string& tmpName,
                       const std::string& key,
                       const std::string& value)
{
    namespace fs = std::filesystem;
    const fs::path p = fs::temp_directory_path() / tmpName;
    {
        std::ofstream ofs(p);
        ofs << key << " = " << value << "\n";
    }
    Config cfg;
    cfg.load(p.string());
    fs::remove(p);
    return cfg;
}

}  // namespace

// ── Constructor: solver.dt ────────────────────────────────────────────────────

TEST(NavierStokesSolverParams, Constructor_DtZero_Throws)
{
    Config cfg = singleKeyConfig("p_dt_zero.cfg", "solver.dt", "0.0");
    EXPECT_THROW(NavierStokesSolver solver(cfg), std::invalid_argument);
}

TEST(NavierStokesSolverParams, Constructor_DtNegative_Throws)
{
    Config cfg = singleKeyConfig("p_dt_neg.cfg", "solver.dt", "-0.001");
    EXPECT_THROW(NavierStokesSolver solver(cfg), std::invalid_argument);
}

// ── Constructor: solver.rho ───────────────────────────────────────────────────

TEST(NavierStokesSolverParams, Constructor_RhoZero_Throws)
{
    Config cfg = singleKeyConfig("p_rho_zero.cfg", "solver.rho", "0.0");
    EXPECT_THROW(NavierStokesSolver solver(cfg), std::invalid_argument);
}

TEST(NavierStokesSolverParams, Constructor_RhoNegative_Throws)
{
    Config cfg = singleKeyConfig("p_rho_neg.cfg", "solver.rho", "-1.0");
    EXPECT_THROW(NavierStokesSolver solver(cfg), std::invalid_argument);
}

// ── Constructor: solver.nu ────────────────────────────────────────────────────

TEST(NavierStokesSolverParams, Constructor_NuNegative_Throws)
{
    Config cfg = singleKeyConfig("p_nu_neg.cfg", "solver.nu", "-0.001");
    EXPECT_THROW(NavierStokesSolver solver(cfg), std::invalid_argument);
}

// nu = 0 is physically meaningful (inviscid limit) and must NOT throw.
TEST(NavierStokesSolverParams, Constructor_NuZero_DoesNotThrow)
{
    Config cfg = singleKeyConfig("p_nu_zero.cfg", "solver.nu", "0.0");
    EXPECT_NO_THROW(NavierStokesSolver solver(cfg));
}

// ── Constructor: solver.tolerance ────────────────────────────────────────────

TEST(NavierStokesSolverParams, Constructor_ToleranceZero_Throws)
{
    Config cfg = singleKeyConfig("p_tol_zero.cfg", "solver.tolerance", "0.0");
    EXPECT_THROW(NavierStokesSolver solver(cfg), std::invalid_argument);
}

TEST(NavierStokesSolverParams, Constructor_ToleranceNegative_Throws)
{
    Config cfg = singleKeyConfig("p_tol_neg.cfg", "solver.tolerance", "-1e-6");
    EXPECT_THROW(NavierStokesSolver solver(cfg), std::invalid_argument);
}

// ── Constructor: solver.alpha_u ───────────────────────────────────────────────

TEST(NavierStokesSolverParams, Constructor_AlphaUZero_Throws)
{
    Config cfg = singleKeyConfig("p_au_zero.cfg", "solver.alpha_u", "0.0");
    EXPECT_THROW(NavierStokesSolver solver(cfg), std::invalid_argument);
}

TEST(NavierStokesSolverParams, Constructor_AlphaUAboveOne_Throws)
{
    Config cfg = singleKeyConfig("p_au_big.cfg", "solver.alpha_u", "1.1");
    EXPECT_THROW(NavierStokesSolver solver(cfg), std::invalid_argument);
}

// alpha_u = 1.0 is the boundary value — must be valid (no under-relaxation).
TEST(NavierStokesSolverParams, Constructor_AlphaUOne_DoesNotThrow)
{
    Config cfg = singleKeyConfig("p_au_one.cfg", "solver.alpha_u", "1.0");
    EXPECT_NO_THROW(NavierStokesSolver solver(cfg));
}

// ── Constructor: solver.alpha_p ───────────────────────────────────────────────

TEST(NavierStokesSolverParams, Constructor_AlphaPZero_Throws)
{
    Config cfg = singleKeyConfig("p_ap_zero.cfg", "solver.alpha_p", "0.0");
    EXPECT_THROW(NavierStokesSolver solver(cfg), std::invalid_argument);
}

TEST(NavierStokesSolverParams, Constructor_AlphaPAboveOne_Throws)
{
    Config cfg = singleKeyConfig("p_ap_big.cfg", "solver.alpha_p", "1.5");
    EXPECT_THROW(NavierStokesSolver solver(cfg), std::invalid_argument);
}

// alpha_p = 1.0 is valid.
TEST(NavierStokesSolverParams, Constructor_AlphaPOne_DoesNotThrow)
{
    Config cfg = singleKeyConfig("p_ap_one.cfg", "solver.alpha_p", "1.0");
    EXPECT_NO_THROW(NavierStokesSolver solver(cfg));
}

// ── Constructor: output.vtk_interval ─────────────────────────────────────────

TEST(NavierStokesSolverParams, Constructor_VtkIntervalZero_Throws)
{
    Config cfg = singleKeyConfig("p_vtki_zero.cfg", "output.vtk_interval", "0");
    EXPECT_THROW(NavierStokesSolver solver(cfg), std::invalid_argument);
}

TEST(NavierStokesSolverParams, Constructor_VtkIntervalNegative_Throws)
{
    Config cfg = singleKeyConfig("p_vtki_neg.cfg", "output.vtk_interval", "-5");
    EXPECT_THROW(NavierStokesSolver solver(cfg), std::invalid_argument);
}

// ── run(): maxIter < 0 ────────────────────────────────────────────────────────

TEST(NavierStokesSolverParams, Run_MaxIterNegative_Throws)
{
    Config cfg;  // all defaults — valid
    NavierStokesSolver solver(cfg);
    solver.initialize();
    // Throw happens before any filesystem operations, so no cleanup needed.
    EXPECT_THROW(solver.run(-1), std::invalid_argument);
}

// ── step(): dt <= 0 ───────────────────────────────────────────────────────────

TEST(NavierStokesSolverParams, Step_DtZeroAfterInit_Throws)
{
    Config cfg;
    NavierStokesSolver solver(cfg);
    solver.initialize();
    EXPECT_THROW(solver.step(0.0), std::invalid_argument);
}

TEST(NavierStokesSolverParams, Step_DtNegativeAfterInit_Throws)
{
    Config cfg;
    NavierStokesSolver solver(cfg);
    solver.initialize();
    EXPECT_THROW(solver.step(-0.01), std::invalid_argument);
}
