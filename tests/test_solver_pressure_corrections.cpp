#include "solver/NavierStokesSolver.hpp"
#include "utils/Config.hpp"

#include <filesystem>
#include <fstream>
#include <gtest/gtest.h>
#include <stdexcept>
#include <string>

// Tests for solver.pressure_corrections_per_step config key.

namespace
{

/// Writes a minimal 4x4 config with an optional extra key, loads it.
Config minimalConfig(const std::string& tmpName,
                     const std::string& extraKey   = "",
                     const std::string& extraValue = "")
{
    namespace fs = std::filesystem;
    const fs::path p = fs::temp_directory_path() / tmpName;
    {
        std::ofstream ofs(p);
        ofs << "mesh.Nx = 4\n"
               "mesh.Ny = 4\n"
               "mesh.Lx = 1.0\n"
               "mesh.Ly = 1.0\n"
               "solver.dt = 0.001\n"
               "solver.nu = 0.01\n"
               "bc.top.type = WALL\n"
               "bc.top.value = 1.0\n"
               "bc.bottom.type = WALL\n"
               "bc.left.type = WALL\n"
               "bc.right.type = WALL\n"
               "output.vtk_interval = 100\n";
        if (!extraKey.empty())
            ofs << extraKey << " = " << extraValue << "\n";
    }
    Config cfg;
    cfg.load(p.string());
    fs::remove(p);
    return cfg;
}

}  // namespace

// ── Parameter validation ───────────────────────────────────────────────────────

TEST(NavierStokesSolverPressureCorrections, ZeroCorrections_Throws)
{
    // pressure_corrections_per_step = 0 is invalid (must be >= 1).
    Config cfg = minimalConfig("pcp_zero.cfg",
                               "solver.pressure_corrections_per_step", "0");
    EXPECT_THROW(NavierStokesSolver solver(cfg), std::invalid_argument);
}

TEST(NavierStokesSolverPressureCorrections, NegativeCorrections_Throws)
{
    Config cfg = minimalConfig("pcp_neg.cfg",
                               "solver.pressure_corrections_per_step", "-1");
    EXPECT_THROW(NavierStokesSolver solver(cfg), std::invalid_argument);
}

TEST(NavierStokesSolverPressureCorrections, OneCorrection_IsDefault_Constructs)
{
    // Key absent → default is 1; must not throw.
    Config cfg = minimalConfig("pcp_default.cfg");
    EXPECT_NO_THROW(NavierStokesSolver solver(cfg));
}

TEST(NavierStokesSolverPressureCorrections, TwoCorrections_Constructs)
{
    Config cfg = minimalConfig("pcp_two.cfg",
                               "solver.pressure_corrections_per_step", "2");
    EXPECT_NO_THROW(NavierStokesSolver solver(cfg));
}

// ── Behaviour: 2 corrections per step must not worsen continuity residual ─────

TEST(NavierStokesSolverPressureCorrections, TwoCorrections_DoNotWorsenContinuityResidual)
{
    // Run two identical solvers for 3 steps, differing only in
    // pressure_corrections_per_step (1 vs 2).
    // The continuity residual with 2 corrections must be <= that with 1.
    // This uses a lid-driven cavity (top WALL u=1) on a 4x4 mesh as the
    // deterministic non-trivial test case.
    Config cfg1 = minimalConfig("pcp_cmp1.cfg",
                                "solver.pressure_corrections_per_step", "1");
    Config cfg2 = minimalConfig("pcp_cmp2.cfg",
                                "solver.pressure_corrections_per_step", "2");

    NavierStokesSolver solver1(cfg1);
    NavierStokesSolver solver2(cfg2);

    solver1.initialize();
    solver2.initialize();

    constexpr int kSteps = 3;
    constexpr double kDt = 0.001;
    for (int s = 0; s < kSteps; ++s)
    {
        solver1.step(kDt);
        solver2.step(kDt);
    }

    const double cont1 = solver1.continuityResidual();
    const double cont2 = solver2.continuityResidual();

    // Both must be finite.
    EXPECT_TRUE(std::isfinite(cont1)) << "1-correction continuity residual is non-finite";
    EXPECT_TRUE(std::isfinite(cont2)) << "2-correction continuity residual is non-finite";

    // 2 corrections per step must not produce a higher continuity residual.
    // Allow a small floating-point tolerance.
    EXPECT_LE(cont2, cont1 * (1.0 + 1e-6))
        << "2 corrections per step worsened continuity residual: "
        << cont2 << " > " << cont1;
}

// ── Behaviour: residuals remain finite with multiple corrections ───────────────

TEST(NavierStokesSolverPressureCorrections, FiveCorrections_ResidualsFinite)
{
    Config cfg = minimalConfig("pcp_five.cfg",
                               "solver.pressure_corrections_per_step", "5");
    NavierStokesSolver solver(cfg);
    solver.initialize();
    EXPECT_NO_THROW(solver.step(0.001));
    EXPECT_TRUE(std::isfinite(solver.velocityResidual()));
    EXPECT_TRUE(std::isfinite(solver.continuityResidual()));
    EXPECT_TRUE(std::isfinite(solver.pressureResidual()));
}
