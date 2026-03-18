#include "solver/NavierStokesSolver.hpp"
#include "utils/Config.hpp"

#include <filesystem>
#include <fstream>
#include <gtest/gtest.h>
#include <stdexcept>
#include <string>

// Tests for convection scheme selection and CFL safety clamp.

// ── Helpers ───────────────────────────────────────────────────────────────────

namespace
{

/// Writes a single key=value config, loads it, deletes the temp file.
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

// ── Convection scheme parsing ─────────────────────────────────────────────────

TEST(NavierStokesSolverConvection, UnknownScheme_ThrowsInvalidArgument)
{
    Config cfg = singleKeyConfig("conv_bad.cfg",
                                 "solver.convection_scheme", "bogus");
    EXPECT_THROW(NavierStokesSolver solver(cfg), std::invalid_argument);
}

TEST(NavierStokesSolverConvection, CentralScheme_Constructs)
{
    Config cfg = minimalConfig("conv_central.cfg",
                               "solver.convection_scheme", "central");
    EXPECT_NO_THROW(NavierStokesSolver solver(cfg));
}

TEST(NavierStokesSolverConvection, UpwindScheme_Constructs)
{
    Config cfg = minimalConfig("conv_upwind.cfg",
                               "solver.convection_scheme", "upwind");
    EXPECT_NO_THROW(NavierStokesSolver solver(cfg));
}

TEST(NavierStokesSolverConvection, CentralScheme_CaseInsensitive)
{
    Config cfg = minimalConfig("conv_central_upper.cfg",
                               "solver.convection_scheme", "CENTRAL");
    EXPECT_NO_THROW(NavierStokesSolver solver(cfg));
}

TEST(NavierStokesSolverConvection, DefaultScheme_IsUpwind_StepProducesFiniteResult)
{
    // No solver.convection_scheme key → default is upwind.
    Config cfg = minimalConfig("conv_default.cfg");
    NavierStokesSolver solver(cfg);
    solver.initialize();
    EXPECT_NO_THROW(solver.step(0.001));
    EXPECT_TRUE(std::isfinite(solver.velocityResidual()));
    EXPECT_TRUE(std::isfinite(solver.continuityResidual()));
}

// ── CFL parameters ────────────────────────────────────────────────────────────

TEST(NavierStokesSolverConvection, MaxCflConv_Zero_Throws)
{
    Config cfg = singleKeyConfig("cfl_conv_zero.cfg", "solver.max_cfl_conv", "0.0");
    EXPECT_THROW(NavierStokesSolver solver(cfg), std::invalid_argument);
}

TEST(NavierStokesSolverConvection, MaxCflConv_Negative_Throws)
{
    Config cfg = singleKeyConfig("cfl_conv_neg.cfg", "solver.max_cfl_conv", "-1.0");
    EXPECT_THROW(NavierStokesSolver solver(cfg), std::invalid_argument);
}

TEST(NavierStokesSolverConvection, MaxCflDiff_Zero_Throws)
{
    Config cfg = singleKeyConfig("cfl_diff_zero.cfg", "solver.max_cfl_diff", "0.0");
    EXPECT_THROW(NavierStokesSolver solver(cfg), std::invalid_argument);
}

// ── dtEffective accessor ──────────────────────────────────────────────────────

TEST(NavierStokesSolverConvection, DtEffective_NeverExceedsUserDt)
{
    // After one step, dtEffective() must be <= the configured dt.
    Config cfg = minimalConfig("cfl_dteff.cfg");
    NavierStokesSolver solver(cfg);
    solver.initialize();
    solver.step(0.001);
    EXPECT_LE(solver.dtEffective(), 0.001 * (1.0 + 1e-9));
    EXPECT_GT(solver.dtEffective(), 0.0);
}

TEST(NavierStokesSolverConvection, UpwindScheme_StepProducesFiniteResult)
{
    Config cfg = minimalConfig("conv_upwind_step.cfg",
                               "solver.convection_scheme", "upwind");
    NavierStokesSolver solver(cfg);
    solver.initialize();
    EXPECT_NO_THROW(solver.step(0.001));
    EXPECT_TRUE(std::isfinite(solver.velocityResidual()));
}

TEST(NavierStokesSolverConvection, CentralScheme_StepProducesFiniteResult)
{
    Config cfg = minimalConfig("conv_central_step.cfg",
                               "solver.convection_scheme", "central");
    NavierStokesSolver solver(cfg);
    solver.initialize();
    EXPECT_NO_THROW(solver.step(0.001));
    EXPECT_TRUE(std::isfinite(solver.velocityResidual()));
}
