#include "core/BoundaryCondition.hpp"
#include "core/Field.hpp"
#include "core/Mesh.hpp"
#include "solver/Discretization.hpp"

#include <Eigen/Dense>
#include <gtest/gtest.h>
#include <stdexcept>
#include <vector>

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

// ── Boundary-aware overloads ──────────────────────────────────────────────────

namespace
{
/// 8×4 channel on [0, 2] × [0, 1]: dx = dy = 0.25.
Mesh makeChannelMesh()
{
    Mesh m;
    m.load(8, 4, 2.0, 1.0);
    return m;
}

}  // namespace

TEST(Discretization, GradientPressureFaces_LinearWithFixedOutlet_ExactInOutletColumn)
{
    const Mesh mesh = makeChannelMesh();
    BoundaryCondition bc;
    bc.addPatch("right", {BoundaryType::OUTLET, Eigen::Vector2d::Zero()});
    const std::vector<PressureBoundaryFace> faces = bc.buildPressureFaces(mesh);

    // p = 2 − x: zero at the outlet face, gradient −1.
    Field<double> p(mesh, 0.0);
    for (int i = 0; i < mesh.Nx(); ++i)
        for (int j = 0; j < mesh.Ny(); ++j)
            p(i, j) = 2.0 - mesh.getCellCenter(i, j).x();

    const Field<Eigen::Vector2d> grad = Discretization::gradient(p, mesh, faces);

    // Interior columns and the outlet column (fixed face value) are exact;
    // the left column has a zero-gradient face and is not.
    for (int i = 1; i < mesh.Nx(); ++i)
        for (int j = 0; j < mesh.Ny(); ++j)
        {
            EXPECT_NEAR(grad(i, j).x(), -1.0, 1e-12) << "i=" << i << " j=" << j;
            EXPECT_NEAR(grad(i, j).y(),  0.0, 1e-12) << "i=" << i << " j=" << j;
        }
}

TEST(Discretization, GradientVelocityFaces_LinearShearBetweenWalls_ExactAtWallRows)
{
    const Mesh mesh = makeChannelMesh();
    BoundaryCondition bc;
    bc.addPatch("bottom", {BoundaryType::WALL, Eigen::Vector2d(0.0, 0.0)});
    bc.addPatch("top",    {BoundaryType::WALL, Eigen::Vector2d(1.0, 0.0)});
    const std::vector<VelocityBoundaryFace> faces = bc.buildVelocityFaces(mesh);

    // Couette profile u_x = y matches both wall values: ∂u_x/∂y = 1 everywhere.
    Field<double> ux(mesh, 0.0);
    for (int i = 0; i < mesh.Nx(); ++i)
        for (int j = 0; j < mesh.Ny(); ++j)
            ux(i, j) = mesh.getCellCenter(i, j).y();

    const Field<Eigen::Vector2d> grad = Discretization::gradient(ux, mesh, faces, 0);

    for (int i = 0; i < mesh.Nx(); ++i)
    {
        EXPECT_NEAR(grad(i, 0).y(),              1.0, 1e-12) << "bottom row i=" << i;
        EXPECT_NEAR(grad(i, mesh.Ny() - 1).y(),  1.0, 1e-12) << "top row i=" << i;
    }
}

TEST(Discretization, GradientVelocityFaces_BadAxis_Throws)
{
    const Mesh mesh = makeChannelMesh();
    const BoundaryCondition bc;
    const std::vector<VelocityBoundaryFace> faces = bc.buildVelocityFaces(mesh);
    const Field<double> ux(mesh, 0.0);

    EXPECT_THROW(static_cast<void>(Discretization::gradient(ux, mesh, faces, 2)),
                 std::invalid_argument);
}

TEST(Discretization, RhieChowFaces_InletIntoStillFluid_ExactInflowInFirstColumn)
{
    const Mesh mesh = makeChannelMesh();
    BoundaryCondition bc;
    bc.addPatch("left",   {BoundaryType::INLET,  Eigen::Vector2d(1.0, 0.0)});
    bc.addPatch("right",  {BoundaryType::OUTLET, Eigen::Vector2d::Zero()});
    bc.addPatch("bottom", {BoundaryType::WALL,   Eigen::Vector2d::Zero()});
    bc.addPatch("top",    {BoundaryType::WALL,   Eigen::Vector2d::Zero()});
    const std::vector<VelocityBoundaryFace> vf = bc.buildVelocityFaces(mesh);
    const std::vector<PressureBoundaryFace> pf = bc.buildPressureFaces(mesh);

    const Field<Eigen::Vector2d> vel(mesh, Eigen::Vector2d::Zero());
    const Field<double>          p(mesh, 0.0);

    const Field<double> div =
        Discretization::divergenceRhieChow(vel, p, mesh, 0.01, 1.0, vf, pf);

    // Inlet face carries U·A into each first-column cell: div = −U·dy/(dx·dy) = −1/dx.
    for (int j = 0; j < mesh.Ny(); ++j)
        EXPECT_NEAR(div(0, j), -4.0, 1e-12) << "j=" << j;
    for (int i = 1; i < mesh.Nx(); ++i)
        for (int j = 0; j < mesh.Ny(); ++j)
            EXPECT_NEAR(div(i, j), 0.0, 1e-12) << "i=" << i << " j=" << j;
}

TEST(Discretization, RhieChowFaces_LinearPressureToOutlet_NoSpuriousFlux)
{
    const Mesh mesh = makeChannelMesh();
    BoundaryCondition bc;
    bc.addPatch("right", {BoundaryType::OUTLET, Eigen::Vector2d::Zero()});
    const std::vector<VelocityBoundaryFace> vf = bc.buildVelocityFaces(mesh);
    const std::vector<PressureBoundaryFace> pf = bc.buildPressureFaces(mesh);

    // A linear pressure has equal compact and cell gradients, so the outlet
    // Rhie-Chow correction must vanish: still fluid stays divergence-free there.
    Field<double> p(mesh, 0.0);
    for (int i = 0; i < mesh.Nx(); ++i)
        for (int j = 0; j < mesh.Ny(); ++j)
            p(i, j) = 2.0 - mesh.getCellCenter(i, j).x();
    const Field<Eigen::Vector2d> vel(mesh, Eigen::Vector2d::Zero());

    const Field<double> div =
        Discretization::divergenceRhieChow(vel, p, mesh, 0.01, 1.0, vf, pf);

    for (int j = 0; j < mesh.Ny(); ++j)
        EXPECT_NEAR(div(mesh.Nx() - 1, j), 0.0, 1e-12) << "outlet j=" << j;
}

TEST(Discretization, RhieChowFaces_MismatchedFaceLists_Throws)
{
    const Mesh mesh = makeChannelMesh();
    const BoundaryCondition bc;
    const std::vector<VelocityBoundaryFace> vf = bc.buildVelocityFaces(mesh);
    std::vector<PressureBoundaryFace>       pf = bc.buildPressureFaces(mesh);
    pf.pop_back();

    const Field<Eigen::Vector2d> vel(mesh, Eigen::Vector2d::Zero());
    const Field<double>          p(mesh, 0.0);

    EXPECT_THROW(static_cast<void>(
                     Discretization::divergenceRhieChow(vel, p, mesh, 0.01, 1.0, vf, pf)),
                 std::invalid_argument);
}
