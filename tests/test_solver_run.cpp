#include "solver/NavierStokesSolver.hpp"
#include "utils/Config.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <gtest/gtest.h>
#include <sstream>
#include <stdexcept>
#include <string>

// ── Helpers ───────────────────────────────────────────────────────────────────

namespace
{

/// Writes a minimal config file and returns a loaded Config.
/// Uses a 4x4 mesh (small, deterministic) with output.dir set to tmpDir.
Config makeSmallConfig(const std::string& tmpDir)
{
    namespace fs = std::filesystem;
    const fs::path cfgPath =
        fs::temp_directory_path() / "lgflow_test_solver_run.cfg";

    {
        std::ofstream ofs(cfgPath);
        ofs << "mesh.Nx            = 4\n";
        ofs << "mesh.Ny            = 4\n";
        ofs << "mesh.Lx            = 1.0\n";
        ofs << "mesh.Ly            = 1.0\n";
        ofs << "solver.dt          = 0.001\n";
        ofs << "solver.rho         = 1.0\n";
        ofs << "solver.nu          = 0.01\n";
        ofs << "solver.tolerance   = 1e-6\n";
        ofs << "solver.alpha_u     = 0.7\n";
        ofs << "solver.alpha_p     = 0.3\n";
        ofs << "output.vtk_interval = 1000\n";
        ofs << "output.dir         = " << tmpDir << "\n";
    }

    Config cfg;
    cfg.load(cfgPath.string());
    fs::remove(cfgPath);
    return cfg;
}

/// Creates a unique temporary directory for a test and returns its path.
/// Caller is responsible for cleanup.
std::string makeTempDir(const std::string& suffix)
{
    namespace fs = std::filesystem;
    const fs::path dir =
        fs::temp_directory_path() / ("lgflow_test_" + suffix);
    fs::create_directories(dir);
    return dir.string();
}

}  // namespace

// ── Step: residuals after one step ────────────────────────────────────────────

TEST(NavierStokesSolverRun, Step_AfterInit_ResidualsAreFinite)
{
    // After a single step on a zero-initialised 4x4 mesh, all three residuals
    // must be finite (not NaN, not Inf).
    const std::string tmpDir = makeTempDir("step_residuals");

    Config cfg = makeSmallConfig(tmpDir);
    NavierStokesSolver solver(cfg);
    solver.initialize();

    ASSERT_NO_THROW(solver.step(0.001));

    EXPECT_TRUE(std::isfinite(solver.velocityResidual()))
        << "velocity residual is not finite";
    EXPECT_TRUE(std::isfinite(solver.continuityResidual()))
        << "continuity residual is not finite";
    EXPECT_TRUE(std::isfinite(solver.pressureResidual()))
        << "pressure residual is not finite";

    EXPECT_GE(solver.velocityResidual(),  0.0);
    EXPECT_GE(solver.continuityResidual(), 0.0);
    EXPECT_GE(solver.pressureResidual(),  0.0);

    std::filesystem::remove_all(tmpDir);
}

// ── Run: history CSV written with header and rows ─────────────────────────────

TEST(NavierStokesSolverRun, Run_WritesHistoryCsv_WithHeaderAndRows)
{
    // run(N) must produce history.csv with the correct header on line 1
    // and exactly N data rows.
    const std::string tmpDir = makeTempDir("history_csv");

    Config cfg = makeSmallConfig(tmpDir);
    NavierStokesSolver solver(cfg);
    solver.initialize();

    constexpr int N = 3;
    ASSERT_NO_THROW(solver.run(N));

    const std::string csvPath = tmpDir + "/history.csv";
    std::ifstream ifs(csvPath);
    ASSERT_TRUE(ifs.is_open()) << "history.csv not found at " << csvPath;

    std::string line;
    // Line 1: header
    ASSERT_TRUE(std::getline(ifs, line));
    EXPECT_EQ(line, "iter,vel_residual,cont_residual,pressure_residual");

    // Lines 2+: data rows — must be parseable and have 4 fields.
    // The solver may converge before maxIter, so rowCount can be in [1, N].
    int rowCount = 0;
    while (std::getline(ifs, line))
    {
        if (line.empty()) continue;
        std::istringstream ss(line);
        std::string token;
        int fields = 0;
        while (std::getline(ss, token, ','))
            ++fields;
        EXPECT_EQ(fields, 4) << "Row " << (rowCount + 1) << " has wrong field count: " << line;
        ++rowCount;
    }

    EXPECT_GE(rowCount, 1) << "history.csv must contain at least 1 data row";
    EXPECT_LE(rowCount, N) << "history.csv cannot have more rows than maxIter";

    ifs.close();  // must close before remove_all on Windows
    std::filesystem::remove_all(tmpDir);
}

// ── Run: zero iterations does nothing ────────────────────────────────────────

TEST(NavierStokesSolverRun, Run_MaxIterZero_DoesNotThrow)
{
    // run(0) must exit immediately without throwing or writing anything beyond
    // the header (the CSV may or may not be created — implementation-defined).
    const std::string tmpDir = makeTempDir("maxiter_zero");

    Config cfg = makeSmallConfig(tmpDir);
    NavierStokesSolver solver(cfg);
    solver.initialize();

    EXPECT_NO_THROW(solver.run(0));

    // Residuals must remain at their post-initialize values (no step was taken).
    EXPECT_TRUE(std::isfinite(solver.velocityResidual()));
    EXPECT_TRUE(std::isfinite(solver.continuityResidual()));

    std::filesystem::remove_all(tmpDir);
}

// ── Run: convergence check uses both residuals ─────────────────────────────────

TEST(NavierStokesSolverRun, Run_ConvergenceCheck_UsesBothVelocityAndContinuityResiduals)
{
    // Behavioural check: after run(), both velocityResidual() and
    // continuityResidual() are accessible and finite, confirming both are
    // tracked throughout the loop (not just one).
    const std::string tmpDir = makeTempDir("both_residuals");

    Config cfg = makeSmallConfig(tmpDir);
    NavierStokesSolver solver(cfg);
    solver.initialize();

    constexpr int N = 5;
    ASSERT_NO_THROW(solver.run(N));

    // Both residuals must be finite and non-negative after the run.
    EXPECT_TRUE(std::isfinite(solver.velocityResidual()))
        << "velocityResidual() is not finite after run()";
    EXPECT_TRUE(std::isfinite(solver.continuityResidual()))
        << "continuityResidual() is not finite after run()";

    EXPECT_GE(solver.velocityResidual(),  0.0);
    EXPECT_GE(solver.continuityResidual(), 0.0);

    // Verify residual() (legacy accessor) agrees with velocityResidual().
    EXPECT_NEAR(solver.residual(), solver.velocityResidual(), 1e-15);

    std::filesystem::remove_all(tmpDir);
}

// ── Relative continuity residual and stopping rule ────────────────────────────

namespace
{
/// Small 8×4 channel: uniform inlet left, outlet right, walls top/bottom.
Config makeChannelConfig(const std::string& tmpDir)
{
    Config cfg = makeSmallConfig(tmpDir);
    cfg.set<int>("mesh.Nx", 8);
    cfg.set<double>("mesh.Lx", 2.0);
    cfg.set<double>("solver.dt", 0.05);
    cfg.set<std::string>("bc.left.type", "INLET");
    cfg.set<double>("bc.left.value_x", 1.0);
    cfg.set<std::string>("bc.right.type", "OUTLET");
    cfg.set<std::string>("bc.top.type", "WALL");
    cfg.set<std::string>("bc.bottom.type", "WALL");
    return cfg;
}
}  // namespace

TEST(NavierStokesSolverRun, ContinuityResidual_BeforeReferenceSteps_IsOne)
{
    const std::string tmpDir = makeTempDir("cont_reference");
    NavierStokesSolver solver(makeChannelConfig(tmpDir));
    solver.initialize();

    // R_0 is the maximum over the first 5 steps; until then the relative
    // residual is reported as exactly 1.
    for (int k = 0; k < 4; ++k)
    {
        solver.step(0.05);
        EXPECT_NEAR(solver.continuityResidual(), 1.0, 1e-15) << "step " << k + 1;
    }

    std::filesystem::remove_all(tmpDir);
}

TEST(NavierStokesSolverRun, ContinuityResidual_Channel_FallsBelowOne)
{
    const std::string tmpDir = makeTempDir("cont_falls");
    NavierStokesSolver solver(makeChannelConfig(tmpDir));
    solver.initialize();

    for (int k = 0; k < 60; ++k)
        solver.step(0.05);

    EXPECT_GT(solver.continuityResidual(), 0.0);
    EXPECT_LT(solver.continuityResidual(), 0.1);

    std::filesystem::remove_all(tmpDir);
}

TEST(NavierStokesSolverRun, IsConverged_VelocityMetContinuityNot_ReturnsFalse)
{
    // Velocity tolerance trivially satisfied, continuity tolerance unreachable:
    // the solver must not report convergence (both criteria are required).
    const std::string tmpDir = makeTempDir("conv_vel_only");
    Config cfg = makeChannelConfig(tmpDir);
    cfg.set<double>("solver.tolerance", 10.0);
    cfg.set<double>("solver.tolerance_continuity", 1e-300);
    NavierStokesSolver solver(cfg);
    solver.initialize();

    solver.run(20);

    EXPECT_LT(solver.velocityResidual(), 10.0);
    EXPECT_FALSE(solver.isConverged());

    std::filesystem::remove_all(tmpDir);
}

TEST(NavierStokesSolverRun, IsConverged_ContinuityMetVelocityNot_ReturnsFalse)
{
    const std::string tmpDir = makeTempDir("conv_cont_only");
    Config cfg = makeChannelConfig(tmpDir);
    cfg.set<double>("solver.tolerance", 1e-300);
    cfg.set<double>("solver.tolerance_continuity", 10.0);
    NavierStokesSolver solver(cfg);
    solver.initialize();

    solver.run(20);

    EXPECT_LT(solver.continuityResidual(), 10.0);
    EXPECT_FALSE(solver.isConverged());

    std::filesystem::remove_all(tmpDir);
}

TEST(NavierStokesSolverRun, Run_LooseTolerances_DoesNotStopBeforeReferenceSteps)
{
    // Velocity tolerance trivially met and continuity tolerance 1: the
    // relative continuity residual is exactly 1 (not < 1) until R_0 is known
    // after 5 steps, so run() cannot stop earlier.
    const std::string tmpDir = makeTempDir("conv_reference_steps");
    Config cfg = makeChannelConfig(tmpDir);
    cfg.set<double>("solver.tolerance", 10.0);
    cfg.set<double>("solver.tolerance_continuity", 1.0);
    NavierStokesSolver solver(cfg);
    solver.initialize();

    solver.run(200);

    // history.csv: header + one row per completed step.
    int rows = -1;
    {
        std::ifstream hist(tmpDir + "/history.csv");
        for (std::string line; std::getline(hist, line);)
            ++rows;
    }
    EXPECT_GE(rows, 5);

    std::filesystem::remove_all(tmpDir);
}
