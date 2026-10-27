/**
 * @brief [SLOW] Plane Poiseuille channel validation with grid refinement.
 *
 * A uniform inflow U develops between two no-slip walls into the analytical
 * fully developed profile
 *
 *   u(y) = 6·U·(y/H)·(1 − y/H),   u_max = 1.5·U,   Q = U·H
 *
 * (White 2006, "Fluid Mechanics" §6-2). The case is solved on successively
 * refined meshes (Ny doubled per level, fixed cell aspect ratio). Per level:
 *   - the solver must converge (velocity AND relative continuity criteria);
 *   - the flow rate Q(x) = Σ u·dy must equal U·H within validation.max_q_error
 *     at mid-channel, the profile station and the last cell column;
 *   - the Q error at x = 1 H (developing flow) and the RMSE of u(y) at the
 *     profile station, relative to u_max, must both be lower than on the
 *     previous level.
 * Q is summed from cell-centre velocities, which match the conserved face
 * fluxes only to O(dx²·∂²u/∂x²). Near the inlet, where the uniform inflow meets
 * the no-slip walls, that term is large on coarse meshes (0.6 % at x = 1 H for
 * Ny = 8, 4 % in the first cell column), so x = 1 H is required to converge
 * rather than to meet the tolerance on every level.
 * After the levels: the observed order p = log2(E_coarse / E_fine) of the
 * finest pair must lie in [min_order, max_order], the finest RMSE must be
 * below validation.max_rel_rmse and the finest Q error at x = 1 H below
 * validation.max_q_error.
 *
 * Any failed level check stops further refinement. validation.min_levels are
 * always run; further levels run while the predicted time
 * t_next = t_last · (t_last / t_prev) fits the remaining budget, set by the
 * environment variable LGFLOW_VALIDATION_BUDGET_S (default 1800 s). Raising it
 * above the CTest TIMEOUT also needs `ctest --timeout <s>`.
 *
 * Domain length: Lx = entrance_length_factor · L_e, rounded up to whole cells
 * of the coarsest mesh so every level has the same Lx. L_e is the parallel-
 * plate hydrodynamic entrance length of Chen (1973), J. Fluids Eng. 95:153,
 *   L_e / D_h = 0.315 / (1 + 0.0175·Re) + 0.011·Re,  D_h = 2H,  Re = U·D_h/ν,
 * coefficients as cited in later entrance-length studies.
 *
 * Run from FlowCore/ so that cases/ and output/ paths resolve correctly.
 * CTest sets WORKING_DIRECTORY automatically via CMakeLists.txt.
 */

#include "core/Field.hpp"
#include "core/Mesh.hpp"
#include "solver/NavierStokesSolver.hpp"
#include "utils/Config.hpp"

#include <Eigen/Dense>
#include <gtest/gtest.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <format>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

constexpr double DEFAULT_BUDGET_S = 1800.0;

/// Results of one refinement level.
struct LevelResult
{
    int                 nx{0};
    int                 ny{0};
    double              seconds{0.0};
    bool                converged{false};
    double              velResidual{0.0};
    double              contResidual{0.0};
    double              relRmse{0.0};
    std::vector<double> qErrors;  ///< |Q − U·H| / (U·H) per station
};

/// Parallel-plate hydrodynamic entrance length, Chen (1973).
double entranceLength(double U, double H, double nu)
{
    const double dh = 2.0 * H;
    const double re = U * dh / nu;
    return dh * (0.315 / (1.0 + 0.0175 * re) + 0.011 * re);
}

/// Time budget in seconds from LGFLOW_VALIDATION_BUDGET_S, else the default.
double budgetSeconds()
{
    const char* env = std::getenv("LGFLOW_VALIDATION_BUDGET_S");
    if (env == nullptr)
        return DEFAULT_BUDGET_S;
    try
    {
        const double value = std::stod(env);
        if (value > 0.0)
            return value;
    }
    catch (const std::exception&)
    {
    }
    std::cerr << "LGFLOW_VALIDATION_BUDGET_S='" << env
              << "' is not a positive number; using " << DEFAULT_BUDGET_S << " s\n";
    return DEFAULT_BUDGET_S;
}

/// u_x at row j, linearly interpolated in x between the bracketing cell columns.
double uxAt(const Field<Eigen::Vector2d>& vel, const Mesh& mesh, double x, int j)
{
    const int    nx = mesh.Nx();
    const double dx = 2.0 * mesh.getCellCenter(0, 0).x();
    const double s  = x / dx - 0.5;  // fractional column index
    const int    i0 = std::clamp(static_cast<int>(std::floor(s)), 0, nx - 2);
    const double w  = std::clamp(s - static_cast<double>(i0), 0.0, 1.0);
    return (1.0 - w) * vel(i0, j).x() + w * vel(i0 + 1, j).x();
}

/// Flow rate per unit depth through the cross-section at x: Σ_j u_x·dy.
double flowRate(const Field<Eigen::Vector2d>& vel, const Mesh& mesh, double x)
{
    const double dy = 2.0 * mesh.getCellCenter(0, 0).y();
    double q = 0.0;
    for (int j = 0; j < mesh.Ny(); ++j)
        q += uxAt(vel, mesh, x, j) * dy;
    return q;
}

/// RMSE of u_x(y) against the Poiseuille profile at x, divided by u_max = 1.5·U.
double relativeRmse(const Field<Eigen::Vector2d>& vel, const Mesh& mesh,
                    double x, double U, double H)
{
    double sum = 0.0;
    for (int j = 0; j < mesh.Ny(); ++j)
    {
        const double eta   = mesh.getCellCenter(0, j).y() / H;
        const double exact = 6.0 * U * eta * (1.0 - eta);
        const double err   = uxAt(vel, mesh, x, j) - exact;
        sum += err * err;
    }
    return std::sqrt(sum / static_cast<double>(mesh.Ny())) / (1.5 * U);
}

}  // namespace

TEST(PoiseuilleChannel, RefinementStudy_UniformInlet_SecondOrderConvergence)
{
    const std::string cfgPath = "cases/channel_flow/case_validate.cfg";
    if (!std::filesystem::exists(cfgPath))
        GTEST_SKIP() << "Validation case not found — run from FlowCore/";

    Config base;
    base.load(cfgPath);

    const double U  = base.get<double>("bc.left.value_x");
    const double H  = base.get<double>("mesh.Ly");
    const double nu = base.get<double>("solver.nu");

    const double factor       = base.get<double>("validation.entrance_length_factor");
    const double sampleFrac   = base.get<double>("validation.sample_fraction");
    const int    baseNy       = base.get<int>("validation.base_ny");
    const double aspect       = base.get<double>("validation.aspect_ratio");
    const int    minLevels    = base.get<int>("validation.min_levels");
    const double maxQError    = base.get<double>("validation.max_q_error");
    const double minOrder     = base.get<double>("validation.min_order");
    const double maxOrder     = base.get<double>("validation.max_order");
    const double maxRelRmse   = base.get<double>("validation.max_rel_rmse");
    const int    maxIter      = base.get<int>("solver.max_iter");
    const std::string outRoot = base.get<std::string>("output.dir");
    const double budget       = budgetSeconds();

    ASSERT_GE(minLevels, 3) << "an observed order needs at least 3 levels";

    // Same Lx on every level: whole cells of the coarsest mesh.
    const double coarseDx = aspect * H / static_cast<double>(baseNy);
    const double Le       = entranceLength(U, H, nu);
    const double Lx       = std::ceil(factor * Le / coarseDx) * coarseDx;
    const double xSample  = sampleFrac * Lx;

    std::cout << std::format(
        "Poiseuille refinement: L_e = {:.3f} H, Lx = {:.2f} H, profile at x = {:.2f} H, "
        "budget {:.0f} s\n", Le / H, Lx / H, xSample / H, budget);
    std::cout << "level    Nx x  Ny    time[s]  vel_res   cont_res  rel_RMSE   "
                 "Q_err(1H, mid, sample, out)\n";

    std::vector<LevelResult> results;
    double elapsed = 0.0;

    for (int level = 0;; ++level)
    {
        // Beyond the guaranteed levels, refine only if the predicted time fits.
        if (level >= minLevels)
        {
            const double tLast = results[results.size() - 1].seconds;
            const double tPrev = results[results.size() - 2].seconds;
            const double tNext = tLast * (tLast / std::max(tPrev, 1e-3));
            if (elapsed + tNext > budget)
            {
                std::cout << std::format(
                    "Stopping refinement: level {} predicted {:.0f} s, {:.0f} s of budget left\n",
                    level, tNext, budget - elapsed);
                break;
            }
        }

        LevelResult r;
        r.ny = baseNy << level;
        const double dy = H / static_cast<double>(r.ny);
        r.nx = static_cast<int>(std::lround(Lx / (aspect * dy)));

        Config cfg = base;
        cfg.set<int>("mesh.Nx", r.nx);
        cfg.set<int>("mesh.Ny", r.ny);
        cfg.set<double>("mesh.Lx", Lx);
        cfg.set<std::string>("output.dir", outRoot + "/level" + std::to_string(level));

        NavierStokesSolver solver(cfg);
        solver.initialize();

        const auto t0 = std::chrono::steady_clock::now();
        solver.run(maxIter);
        r.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        elapsed += r.seconds;

        r.converged    = solver.isConverged();
        r.velResidual  = solver.velocityResidual();
        r.contResidual = solver.continuityResidual();

        const Field<Eigen::Vector2d>& vel  = solver.velocity();
        const Mesh&                   mesh = vel.mesh();
        const double dx = aspect * dy;
        for (const double x : {1.0 * H, 0.5 * Lx, xSample, Lx - 0.5 * dx})
            r.qErrors.push_back(std::abs(flowRate(vel, mesh, x) - U * H) / (U * H));
        r.relRmse = relativeRmse(vel, mesh, xSample, U, H);

        std::cout << std::format(
            "{:5d} {:5d} x {:3d} {:10.1f}  {:8.2e}  {:8.2e}  {:8.2e}   "
            "{:.1e}, {:.1e}, {:.1e}, {:.1e}\n",
            level, r.nx, r.ny, r.seconds, r.velResidual, r.contResidual, r.relRmse,
            r.qErrors[0], r.qErrors[1], r.qErrors[2], r.qErrors[3]);

        // ── Per-level objectives ─────────────────────────────────────────────
        EXPECT_TRUE(r.converged)
            << "level " << level << " did not converge in " << maxIter << " iterations";
        // Station 0 (x = 1 H) is checked for convergence below, not against the tolerance.
        const char* stations[] = {"x = 1 H", "mid-channel", "profile station", "outlet column"};
        for (std::size_t k = 1; k < r.qErrors.size(); ++k)
            EXPECT_LE(r.qErrors[k], maxQError)
                << "level " << level << ": Q error at " << stations[k];
        if (!results.empty())
        {
            EXPECT_LT(r.qErrors[0], results.back().qErrors[0])
                << "level " << level << ": Q error at x = 1 H did not decrease under refinement";
            EXPECT_LT(r.relRmse, results.back().relRmse)
                << "level " << level << ": error did not decrease under refinement";
        }

        results.push_back(r);

        // Any failed objective fails the test and stops further refinement.
        if (HasFailure())
            return;
    }

    // ── Objectives on the finest levels ──────────────────────────────────────
    ASSERT_GE(results.size(), static_cast<std::size_t>(minLevels));
    const LevelResult& fine   = results[results.size() - 1];
    const LevelResult& coarse = results[results.size() - 2];
    const double order = std::log(coarse.relRmse / fine.relRmse) / std::log(2.0);

    std::cout << std::format("Observed order (finest pair, Ny {} -> {}): p = {:.3f}\n",
                             coarse.ny, fine.ny, order);

    EXPECT_GE(order, minOrder) << "observed order below the expected second order";
    EXPECT_LE(order, maxOrder) << "observed order implausibly high (iteration error floor?)";
    EXPECT_LT(fine.relRmse, maxRelRmse) << "finest-level error too large";
    EXPECT_LE(fine.qErrors[0], maxQError) << "finest-level Q error at x = 1 H too large";
}
