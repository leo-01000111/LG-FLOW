#include "core/Field.hpp"
#include "core/Mesh.hpp"
#include "solver/Discretization.hpp"

#include <Eigen/Dense>
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

// ── Divergence tests ──────────────────────────────────────────────────────────

TEST(Discretization, DivergenceZeroField_IsZeroEverywhere)
{
    Mesh mesh = make4x4Mesh();
    Field<Eigen::Vector2d> u(mesh, Eigen::Vector2d::Zero());
    Field<double> div = Discretization::divergence(u, mesh);
    for (int k = 0; k < div.size(); ++k)
        EXPECT_NEAR(div[k], 0.0, 1e-12);
}

TEST(Discretization, DivergenceLinearField_ApproxTwoOnInteriorCells)
{
    // u = (x, y)  →  ∇·u = 2 everywhere (analytically)
    // Central differencing gives exact result on interior cells for linear u.
    Mesh mesh = make4x4Mesh();
    Field<Eigen::Vector2d> u(mesh, Eigen::Vector2d::Zero());
    for (int i = 0; i < mesh.Nx(); ++i)
        for (int j = 0; j < mesh.Ny(); ++j)
            u(i, j) = mesh.getCellCenter(i, j);

    Field<double> div = Discretization::divergence(u, mesh);

    // Interior cells: i ∈ [1, Nx-2], j ∈ [1, Ny-2]
    for (int i = 1; i < mesh.Nx() - 1; ++i)
        for (int j = 1; j < mesh.Ny() - 1; ++j)
            EXPECT_NEAR(div(i, j), 2.0, 1e-10);
}

// ── Gradient tests ────────────────────────────────────────────────────────────

TEST(Discretization, GradientConstantScalar_IsZeroEverywhere)
{
    // Gauss identity: Σ_f n_f A_f = 0 for any closed CV → ∇(const) = 0
    Mesh mesh = make4x4Mesh();
    Field<double> phi(mesh, 5.0);
    Field<Eigen::Vector2d> grad = Discretization::gradient(phi, mesh);
    for (int k = 0; k < grad.size(); ++k)
    {
        EXPECT_NEAR(grad[k].x(), 0.0, 1e-12);
        EXPECT_NEAR(grad[k].y(), 0.0, 1e-12);
    }
}

TEST(Discretization, GradientLinearScalar_ApproxExactOnInteriorCells)
{
    // φ = 2x − 3y  →  ∇φ = (2, −3) exactly
    // Central differencing is exact for linear fields.
    Mesh mesh = make4x4Mesh();
    Field<double> phi(mesh, 0.0);
    for (int i = 0; i < mesh.Nx(); ++i)
        for (int j = 0; j < mesh.Ny(); ++j)
        {
            const Eigen::Vector2d c = mesh.getCellCenter(i, j);
            phi(i, j) = 2.0 * c.x() - 3.0 * c.y();
        }

    Field<Eigen::Vector2d> grad = Discretization::gradient(phi, mesh);

    for (int i = 1; i < mesh.Nx() - 1; ++i)
        for (int j = 1; j < mesh.Ny() - 1; ++j)
        {
            EXPECT_NEAR(grad(i, j).x(),  2.0, 1e-10);
            EXPECT_NEAR(grad(i, j).y(), -3.0, 1e-10);
        }
}

// ── Laplacian tests ───────────────────────────────────────────────────────────

TEST(Discretization, LaplacianConstantScalar_IsZeroEverywhere)
{
    // ∂φ/∂n = 0 at every face for φ = const → ∇²φ = 0 everywhere
    Mesh mesh = make4x4Mesh();
    Field<double> phi(mesh, 7.0);
    Field<double> lap = Discretization::laplacian(phi, mesh);
    for (int k = 0; k < lap.size(); ++k)
        EXPECT_NEAR(lap[k], 0.0, 1e-12);
}

TEST(Discretization, LaplacianLinearScalar_IsApproxZeroOnInteriorCells)
{
    // ∇²(ax + by) = 0; two-point stencil is exact for linear functions.
    Mesh mesh = make4x4Mesh();
    Field<double> phi(mesh, 0.0);
    for (int i = 0; i < mesh.Nx(); ++i)
        for (int j = 0; j < mesh.Ny(); ++j)
        {
            const Eigen::Vector2d c = mesh.getCellCenter(i, j);
            phi(i, j) = 2.0 * c.x() - 3.0 * c.y();
        }

    Field<double> lap = Discretization::laplacian(phi, mesh);

    for (int i = 1; i < mesh.Nx() - 1; ++i)
        for (int j = 1; j < mesh.Ny() - 1; ++j)
            EXPECT_NEAR(lap(i, j), 0.0, 1e-10);
}

TEST(Discretization, LaplacianQuadraticScalar_ApproxFourOnInteriorCells)
{
    // φ = x² + y²  →  ∇²φ = 4 exactly (second-order stencil is exact here)
    Mesh mesh = make4x4Mesh();
    Field<double> phi(mesh, 0.0);
    for (int i = 0; i < mesh.Nx(); ++i)
        for (int j = 0; j < mesh.Ny(); ++j)
        {
            const Eigen::Vector2d c = mesh.getCellCenter(i, j);
            phi(i, j) = c.squaredNorm();  // x^2 + y^2
        }

    Field<double> lap = Discretization::laplacian(phi, mesh);

    for (int i = 1; i < mesh.Nx() - 1; ++i)
        for (int j = 1; j < mesh.Ny() - 1; ++j)
            EXPECT_NEAR(lap(i, j), 4.0, 1e-10);
}

// ── Mesh mismatch tests ───────────────────────────────────────────────────────

TEST(Discretization, DivergenceMeshMismatch_Throws)
{
    Mesh mesh1; mesh1.load(4, 4, 1.0, 1.0);
    Mesh mesh2; mesh2.load(4, 4, 1.0, 1.0);
    Field<Eigen::Vector2d> u(mesh1, Eigen::Vector2d::Zero());
    EXPECT_THROW(static_cast<void>(Discretization::divergence(u, mesh2)),
                 std::invalid_argument);
}

TEST(Discretization, GradientMeshMismatch_Throws)
{
    Mesh mesh1; mesh1.load(4, 4, 1.0, 1.0);
    Mesh mesh2; mesh2.load(4, 4, 1.0, 1.0);
    Field<double> phi(mesh1, 0.0);
    EXPECT_THROW(static_cast<void>(Discretization::gradient(phi, mesh2)),
                 std::invalid_argument);
}

TEST(Discretization, LaplacianMeshMismatch_Throws)
{
    Mesh mesh1; mesh1.load(4, 4, 1.0, 1.0);
    Mesh mesh2; mesh2.load(4, 4, 1.0, 1.0);
    Field<double> phi(mesh1, 0.0);
    EXPECT_THROW(static_cast<void>(Discretization::laplacian(phi, mesh2)),
                 std::invalid_argument);
}
