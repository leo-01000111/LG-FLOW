#include "core/Field.hpp"
#include "core/Mesh.hpp"
#include "solver/Discretization.hpp"
#include "solver/PressureSolver.hpp"

#include <Eigen/Dense>
#include <cmath>
#include <gtest/gtest.h>
#include <stdexcept>

namespace
{
Mesh make4x4Mesh()
{
    Mesh m;
    m.load(4, 4, 1.0, 1.0);
    return m;
}
}  // namespace

// ── Input validation ──────────────────────────────────────────────────────────

TEST(PressureSolver, Solve_InvalidDt_Throws)
{
    Mesh mesh = make4x4Mesh();
    PressureSolver solver(mesh);
    Field<Eigen::Vector2d> vel(mesh, Eigen::Vector2d::Zero());
    Field<double>          p(mesh, 0.0);

    EXPECT_THROW(static_cast<void>(solver.solve(vel, p, 0.0, 1.0, 0.5)),
                 std::invalid_argument);
    EXPECT_THROW(static_cast<void>(solver.solve(vel, p, -0.01, 1.0, 0.5)),
                 std::invalid_argument);
}

TEST(PressureSolver, Solve_InvalidRho_Throws)
{
    Mesh mesh = make4x4Mesh();
    PressureSolver solver(mesh);
    Field<Eigen::Vector2d> vel(mesh, Eigen::Vector2d::Zero());
    Field<double>          p(mesh, 0.0);

    EXPECT_THROW(static_cast<void>(solver.solve(vel, p, 0.01, 0.0, 0.5)),
                 std::invalid_argument);
    EXPECT_THROW(static_cast<void>(solver.solve(vel, p, 0.01, -1.0, 0.5)),
                 std::invalid_argument);
}

TEST(PressureSolver, Solve_InvalidAlphaP_Throws)
{
    Mesh mesh = make4x4Mesh();
    PressureSolver solver(mesh);
    Field<Eigen::Vector2d> vel(mesh, Eigen::Vector2d::Zero());
    Field<double>          p(mesh, 0.0);

    // alphaP = 0 is invalid (must be strictly positive)
    EXPECT_THROW(static_cast<void>(solver.solve(vel, p, 0.01, 1.0, 0.0)),
                 std::invalid_argument);
    // alphaP > 1 is invalid
    EXPECT_THROW(static_cast<void>(solver.solve(vel, p, 0.01, 1.0, 1.01)),
                 std::invalid_argument);
    // alphaP = 1.0 is valid (upper boundary)
    EXPECT_NO_THROW(static_cast<void>(solver.solve(vel, p, 0.01, 1.0, 1.0)));
}

TEST(PressureSolver, Solve_MeshMismatch_Throws)
{
    Mesh mesh1 = make4x4Mesh();
    Mesh mesh2;
    mesh2.load(4, 4, 1.0, 1.0);  // different object, same geometry

    PressureSolver solver(mesh1);

    // Velocity field on wrong mesh
    Field<Eigen::Vector2d> velWrong(mesh2, Eigen::Vector2d::Zero());
    Field<double>          pGood(mesh1, 0.0);
    EXPECT_THROW(static_cast<void>(solver.solve(velWrong, pGood, 0.01, 1.0, 0.5)),
                 std::invalid_argument);

    // Pressure field on wrong mesh
    Field<Eigen::Vector2d> velGood(mesh1, Eigen::Vector2d::Zero());
    Field<double>          pWrong(mesh2, 0.0);
    EXPECT_THROW(static_cast<void>(solver.solve(velGood, pWrong, 0.01, 1.0, 0.5)),
                 std::invalid_argument);
}

// ── Physics tests ─────────────────────────────────────────────────────────────

TEST(PressureSolver, Solve_ZeroVelocity_LeavesFieldsNearUnchanged_ResidualNearZero)
{
    // u* = 0 everywhere → div(u*) = 0 → b = 0 → p' = 0.
    // Pressure must remain zero; velocity must remain zero; residual must be ~0.
    Mesh mesh = make4x4Mesh();
    PressureSolver solver(mesh);

    Field<Eigen::Vector2d> vel(mesh, Eigen::Vector2d::Zero());
    Field<double>          p(mesh, 0.0);

    const double residual = solver.solve(vel, p, 0.01, 1.0, 0.5);

    EXPECT_LT(residual, 1e-10);

    for (int k = 0; k < mesh.numCells(); ++k)
    {
        EXPECT_NEAR(p[k], 0.0, 1e-12)
            << "Pressure changed at cell " << k << " despite zero divergence";
        EXPECT_NEAR(vel[k].norm(), 0.0, 1e-12)
            << "Velocity changed at cell " << k << " despite zero gradient";
    }
}

TEST(PressureSolver, Solve_NonZeroDivergence_ReturnsFiniteResidual_UpdatesPressure)
{
    // u* = (x, 0) → div(u*) = 1 on interior cells, non-zero on boundary.
    // The pressure correction must be non-trivial and pressure must be updated.
    Mesh mesh = make4x4Mesh();
    PressureSolver solver(mesh);

    Field<Eigen::Vector2d> vel(mesh, Eigen::Vector2d::Zero());
    Field<double>          p(mesh, 0.0);

    // Set u_x = x (linear in x) so that divergence ≈ 1 everywhere
    for (int i = 0; i < mesh.Nx(); ++i)
        for (int j = 0; j < mesh.Ny(); ++j)
        {
            const Eigen::Vector2d c = mesh.getCellCenter(i, j);
            vel(i, j) = Eigen::Vector2d{c.x(), 0.0};
        }

    const double residual = solver.solve(vel, p, 0.01, 1.0, 0.5);

    // Residual must be finite and non-negative
    EXPECT_TRUE(std::isfinite(residual));
    EXPECT_GE(residual, 0.0);

    // At least some cells must have non-zero pressure after correction
    double pNorm = 0.0;
    for (int k = 0; k < mesh.numCells(); ++k)
        pNorm += p[k] * p[k];
    pNorm = std::sqrt(pNorm);

    EXPECT_GT(pNorm, 1e-8) << "Pressure field was not updated despite non-zero divergence";
}

TEST(PressureSolver, Solve_NonZeroDivergence_ReducesDivergenceNorm)
{
    // After one pressure-correction step the divergence of the corrected
    // velocity should be strictly smaller than before the step.
    // This is the primary invariant of the SIMPLE pressure step.
    Mesh mesh = make4x4Mesh();
    PressureSolver solver(mesh);

    Field<Eigen::Vector2d> vel(mesh, Eigen::Vector2d::Zero());
    Field<double>          p(mesh, 0.0);

    for (int i = 0; i < mesh.Nx(); ++i)
        for (int j = 0; j < mesh.Ny(); ++j)
        {
            const Eigen::Vector2d c = mesh.getCellCenter(i, j);
            vel(i, j) = Eigen::Vector2d{c.x(), 0.0};
        }

    const double divNormBefore = Discretization::divergence(vel, mesh).norm();
    ASSERT_GT(divNormBefore, 1e-10) << "Pre-condition: velocity must have non-zero divergence";

    static_cast<void>(solver.solve(vel, p, 0.01, 1.0, 0.5));

    const double divNormAfter = Discretization::divergence(vel, mesh).norm();
    EXPECT_LT(divNormAfter, divNormBefore)
        << "Divergence norm after pressure correction (" << divNormAfter
        << ") is not less than before (" << divNormBefore << ")";
}

TEST(PressureSolver, Solve_StencilConsistentCorrection_FiniteOnLargerMesh)
{
    // Verifies the inline structured-grid gradient in PressureSolver::solve
    // (which replaces the generic Discretization::gradient loop) produces a
    // finite, non-negative solver residual on a 6x6 mesh, exercising boundary,
    // interior, and corner cell paths.
    //
    // Note: divergence reduction is NOT checked here because the reference-cell
    // pin (p'[0]=0) means the Poisson equation does not drive divergence at
    // cell (0,0) to zero; whether the L2-norm decreases depends on the p'
    // profile. Divergence reduction is covered by the 4x4 test above.
    Mesh mesh;
    mesh.load(6, 6, 1.0, 1.0);
    PressureSolver solver(mesh);

    Field<Eigen::Vector2d> vel(mesh, Eigen::Vector2d::Zero());
    Field<double>          p(mesh, 0.0);

    for (int i = 0; i < mesh.Nx(); ++i)
        for (int j = 0; j < mesh.Ny(); ++j)
        {
            const Eigen::Vector2d c = mesh.getCellCenter(i, j);
            vel(i, j) = Eigen::Vector2d{c.x(), 0.0};
        }

    const double res = solver.solve(vel, p, 0.01, 1.0, 0.5);

    // Solver must produce a finite, non-negative residual.
    EXPECT_TRUE(std::isfinite(res));
    EXPECT_GE(res, 0.0);

    // Pressure must be updated (non-trivial correction).
    double pNorm = 0.0;
    for (int k = 0; k < mesh.numCells(); ++k)
        pNorm += p[k] * p[k];
    EXPECT_GT(std::sqrt(pNorm), 1e-8)
        << "Pressure was not updated on 6x6 mesh";

    // Velocity must remain finite everywhere.
    for (int k = 0; k < mesh.numCells(); ++k)
        EXPECT_TRUE(vel[k].allFinite())
            << "Velocity became non-finite at cell " << k;
}
