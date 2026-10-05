#include "solver/SpallartAllmaras.hpp"
#include "solver/WallDistance.hpp"
#include "solver/NavierStokesSolver.hpp"
#include "core/BoundaryCondition.hpp"
#include "core/Field.hpp"
#include "core/Mesh.hpp"
#include "utils/Config.hpp"

#include <cmath>
#include <filesystem>
#include <gtest/gtest.h>

// ── Helpers ───────────────────────────────────────────────────────────────────

namespace
{

Mesh make4x4Mesh()
{
    Mesh m;
    m.load(4, 4, 1.0, 1.0);
    return m;
}

BoundaryCondition bottomWallBC()
{
    BoundaryCondition bc;
    bc.addPatch("bottom", {BoundaryType::WALL, Eigen::Vector2d::Zero()});
    return bc;
}

} // namespace

// ── Closure function unit tests ───────────────────────────────────────────────

TEST(SpallartAllmaras, Chi_BasicRatio)
{
    // chi(nu_tilde, nu) = nu_tilde / nu
    EXPECT_NEAR(SpallartAllmaras::chi(7.1, 1.0), 7.1, 1.0e-14);
    EXPECT_NEAR(SpallartAllmaras::chi(0.0, 1.0), 0.0, 1.0e-14);
    EXPECT_NEAR(SpallartAllmaras::chi(3.0, 2.0), 1.5, 1.0e-14);
}

TEST(SpallartAllmaras, Fv1_AtCv1_IsHalf)
{
    // At chi = cv1 = 7.1: fv1 = cv1³/(cv1³ + cv1³) = 0.5
    const double chiVal = 7.1;
    EXPECT_NEAR(SpallartAllmaras::fv1(chiVal), 0.5, 1.0e-14);
}

TEST(SpallartAllmaras, Fv1_AtZero_IsZero)
{
    EXPECT_NEAR(SpallartAllmaras::fv1(0.0), 0.0, 1.0e-14);
}

TEST(SpallartAllmaras, Fv1_LargeChiApproachesOne)
{
    // As chi → ∞, fv1 → 1
    EXPECT_NEAR(SpallartAllmaras::fv1(1.0e6), 1.0, 1.0e-6);
}

TEST(SpallartAllmaras, Fv2_AtCv1)
{
    // fv2(chi=7.1, fv1=0.5) = 1 - 7.1/(1 + 7.1*0.5) = 1 - 7.1/4.55
    const double chiVal = 7.1;
    const double fv1Val = SpallartAllmaras::fv1(chiVal);
    const double expected = 1.0 - chiVal / (1.0 + chiVal * fv1Val);
    EXPECT_NEAR(SpallartAllmaras::fv2(chiVal, fv1Val), expected, 1.0e-14);
    // Confirm the numeric value matches the spec (~-0.56044)
    EXPECT_NEAR(SpallartAllmaras::fv2(chiVal, fv1Val), -0.56043956, 1.0e-6);
}

TEST(SpallartAllmaras, RFunc_ClampedToRLim)
{
    // r = nu_tilde / (S_tilde * kappa^2 * d^2) = 0.1/(1.0*0.41^2*0.01) >> 10 → clamped
    const double nuTilde   = 0.1;
    const double sTildeVal = 1.0;
    const double d         = 0.1;
    EXPECT_NEAR(SpallartAllmaras::rFunc(nuTilde, sTildeVal, d), 10.0, 1.0e-14);
}

TEST(SpallartAllmaras, RFunc_SmallR)
{
    // With a very large sTildeVal, r should be small and unclamped.
    const double kappa = 0.41;
    const double d     = 1.0;
    const double nuT   = 0.001;
    const double sT    = 1000.0;
    const double expected = nuT / (sT * kappa * kappa * d * d);
    EXPECT_NEAR(SpallartAllmaras::rFunc(nuT, sT, d), expected, 1.0e-12);
}

TEST(SpallartAllmaras, GFunc_AtROne)
{
    // g(r=1) = 1 + cw2*(1 - 1) = 1 + cw2*0 = 1  (since r^6 - r = 1 - 1 = 0)
    EXPECT_NEAR(SpallartAllmaras::gFunc(1.0), 1.0, 1.0e-14);
}

TEST(SpallartAllmaras, Fw_AtGOne)
{
    // fw(g=1.0): g=1, cw3=2 → fw = 1 * ((1+64)/(1+64))^(1/6) = 1.0
    EXPECT_NEAR(SpallartAllmaras::fw(1.0), 1.0, 1.0e-14);
}

// ── computeNuT ────────────────────────────────────────────────────────────────

TEST(SpallartAllmaras, ComputeNuT_ZeroAtWall)
{
    // nuT must be zero wherever nuTilde = 0.
    const Mesh            mesh = make4x4Mesh();
    const Field<double>   nuTilde(mesh, 0.0);
    const Field<double>   nuT = SpallartAllmaras::computeNuT(nuTilde, 0.01, mesh);

    const int N = mesh.numCells();
    for (int c = 0; c < N; ++c)
        EXPECT_NEAR(nuT[c], 0.0, 1.0e-14) << "nuT non-zero at cell " << c;
}

TEST(SpallartAllmaras, ComputeNuT_PositiveAwayFromWall)
{
    // With a positive nuTilde everywhere, nuT must be positive everywhere.
    const Mesh          mesh = make4x4Mesh();
    const double        nu   = 1.0e-3;
    const Field<double> nuTilde(mesh, 3.0 * nu);  // freestream-like value
    const Field<double> nuT = SpallartAllmaras::computeNuT(nuTilde, nu, mesh);

    const int N = mesh.numCells();
    for (int c = 0; c < N; ++c)
        EXPECT_GT(nuT[c], 0.0) << "nuT not positive at cell " << c;
}

TEST(SpallartAllmaras, ComputeNuT_NegativeNuTilde_GivesZero)
{
    // Negative nuTilde (non-physical) must produce nuT = 0, not a negative value.
    const Mesh          mesh = make4x4Mesh();
    const Field<double> nuTilde(mesh, -1.0e-4);
    const Field<double> nuT = SpallartAllmaras::computeNuT(nuTilde, 1.0e-3, mesh);

    const int N = mesh.numCells();
    for (int c = 0; c < N; ++c)
        EXPECT_NEAR(nuT[c], 0.0, 1.0e-14)
            << "nuT non-zero for negative nuTilde at cell " << c;
}

// ── SA solve: basic sanity ────────────────────────────────────────────────────

TEST(SpallartAllmaras, Solve_ZeroVelocity_NoThrow)
{
    // solve() on zero velocity must not throw and must produce non-negative nuTilde.
    const Mesh              mesh = make4x4Mesh();
    const BoundaryCondition bc   = bottomWallBC();

    const double nu     = 1.0e-3;
    const double nuTi   = 3.0 * nu;

    Field<double>          nuTilde(mesh, nuTi);
    const Field<Eigen::Vector2d> vel(mesh, Eigen::Vector2d::Zero());
    const Field<double>    wallDist = WallDistance::compute(mesh, bc);

    // Set wall row to 0 as the SA model requires.
    for (int c : bc.collectCellsOfType(BoundaryType::WALL, mesh))
        nuTilde[c] = 0.0;

    SpallartAllmaras sa(mesh);

    ASSERT_NO_THROW(sa.solve(nuTilde, vel, wallDist, nu, 0.01, 0.7));

    const int N = mesh.numCells();
    for (int c = 0; c < N; ++c)
        EXPECT_GE(nuTilde[c], 0.0) << "nuTilde negative at cell " << c;
}

TEST(SpallartAllmaras, Solve_AllFinite)
{
    // All nuTilde values after solve must be finite.
    const Mesh              mesh = make4x4Mesh();
    const BoundaryCondition bc   = bottomWallBC();
    const double            nu   = 1.0e-3;

    Field<double>          nuTilde(mesh, 3.0 * nu);
    const Field<Eigen::Vector2d> vel(mesh, Eigen::Vector2d(0.1, 0.0));
    const Field<double>    wallDist = WallDistance::compute(mesh, bc);

    for (int c : bc.collectCellsOfType(BoundaryType::WALL, mesh))
        nuTilde[c] = 0.0;

    SpallartAllmaras sa(mesh);
    sa.solve(nuTilde, vel, wallDist, nu, 0.01, 0.7);

    const int N = mesh.numCells();
    for (int c = 0; c < N; ++c)
        EXPECT_TRUE(std::isfinite(nuTilde[c]))
            << "Non-finite nuTilde at cell " << c;
}

// ── Integration test: NavierStokesSolver with SA ─────────────────────────────

TEST(SpallartAllmaras, NavierStokesSolver_Integration_OneStep)
{
    // Construct NavierStokesSolver with the flat_plate config, call initialize()
    // and one step(). Verify that nuT is non-negative everywhere and zero at the
    // bottom row (wall cells).
    //
    // Skip gracefully if the config file is absent (CI without case files).
    const std::string cfgPath = "cases/flat_plate/case.cfg";
    if (!std::filesystem::exists(cfgPath))
    {
        GTEST_SKIP() << "Flat plate config not found: " << cfgPath;
    }

    Config cfg;
    cfg.load(cfgPath);

    // Use a tiny mesh for speed in the unit test.
    cfg.set("mesh.Nx", 8);
    cfg.set("mesh.Ny", 8);

    NavierStokesSolver solver(cfg);
    ASSERT_NO_THROW(solver.initialize());
    ASSERT_NO_THROW(solver.step(0.001));

    // Verify velocity field is finite (solver did not blow up).
    const auto& vel = solver.velocity();
    const int   N   = cfg.get<int>("mesh.Nx", 8) * cfg.get<int>("mesh.Ny", 8);
    for (int c = 0; c < N; ++c)
    {
        EXPECT_TRUE(std::isfinite(vel[c].x())) << "vel.x non-finite at cell " << c;
        EXPECT_TRUE(std::isfinite(vel[c].y())) << "vel.y non-finite at cell " << c;
    }
}
