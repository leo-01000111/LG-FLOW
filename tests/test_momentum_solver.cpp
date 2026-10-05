#include "solver/MomentumSolver.hpp"
#include "solver/Discretization.hpp"
#include "core/BoundaryCondition.hpp"
#include "core/Mesh.hpp"
#include "core/Field.hpp"

#include <cmath>
#include <gtest/gtest.h>
#include <stdexcept>
#include <utility>
#include <vector>

// ── Helpers ───────────────────────────────────────────────────────────────────

namespace
{

/// Builds a loaded 4×4 unit-square mesh.
Mesh make4x4Mesh()
{
    Mesh m;
    m.load(4, 4, 1.0, 1.0);
    return m;
}

/// Builds a zero velocity field on the given mesh.
Field<Eigen::Vector2d> zeroVelocity(const Mesh& mesh)
{
    return Field<Eigen::Vector2d>(mesh, Eigen::Vector2d::Zero());
}

/// Builds a zero pressure field on the given mesh.
Field<double> zeroPressure(const Mesh& mesh)
{
    return Field<double>(mesh, 0.0);
}

}  // namespace

// ── MomentumSolver: basic sanity ─────────────────────────────────────────────

TEST(MomentumSolver, Solve_ZeroInitial_NoThrow)
{
    // solve() on an all-zero state must not throw and must return a finite uStar.
    Mesh mesh = make4x4Mesh();
    MomentumSolver ms(mesh);

    Field<Eigen::Vector2d> uStar   = zeroVelocity(mesh);
    const auto             uOld    = zeroVelocity(mesh);
    const auto             pressure = zeroPressure(mesh);

    ASSERT_NO_THROW(ms.solve(uStar, uOld, pressure,
                             /*dt=*/0.01, /*rho=*/1.0, /*nu=*/0.01,
                             ConvectionScheme::UPWIND));

    const int N = mesh.numCells();
    for (int c = 0; c < N; ++c)
    {
        EXPECT_TRUE(std::isfinite(uStar[c].x()))
            << "uStar.x is not finite at cell " << c;
        EXPECT_TRUE(std::isfinite(uStar[c].y()))
            << "uStar.y is not finite at cell " << c;
    }
}

TEST(MomentumSolver, Solve_ZeroInitialCentral_NoThrow)
{
    // Same as above but with CENTRAL convection scheme.
    Mesh mesh = make4x4Mesh();
    MomentumSolver ms(mesh);

    Field<Eigen::Vector2d> uStar    = zeroVelocity(mesh);
    const auto             uOld     = zeroVelocity(mesh);
    const auto             pressure = zeroPressure(mesh);

    ASSERT_NO_THROW(ms.solve(uStar, uOld, pressure,
                             0.01, 1.0, 0.01, ConvectionScheme::CENTRAL));

    const int N = mesh.numCells();
    for (int c = 0; c < N; ++c)
    {
        EXPECT_TRUE(std::isfinite(uStar[c].x()));
        EXPECT_TRUE(std::isfinite(uStar[c].y()));
    }
}

// ── MomentumSolver: nu = 0, no pressure → u* = u^k ──────────────────────────

TEST(MomentumSolver, Solve_NuZeroZeroPressure_UStarEqualsUOld)
{
    // When nu = 0 and p = 0, the implicit system reduces to:
    //   (V/dt) u*_P = (V/dt) u^k_P  (diagonal, no off-diagonals)
    // → u* = u^k exactly.
    // Also requires zero convection: use uniform velocity so (u·∇)u = 0.
    Mesh mesh = make4x4Mesh();
    MomentumSolver ms(mesh);

    // Uniform velocity (1, 0) — has zero (u·∇)u.
    Field<Eigen::Vector2d> uOld(mesh, Eigen::Vector2d(1.0, 0.0));
    Field<Eigen::Vector2d> uStar   = zeroVelocity(mesh);
    const auto             pressure = zeroPressure(mesh);

    ms.solve(uStar, uOld, pressure, 0.01, 1.0, /*nu=*/0.0, ConvectionScheme::UPWIND);

    const int N = mesh.numCells();
    for (int c = 0; c < N; ++c)
    {
        EXPECT_NEAR(uStar[c].x(), 1.0, 1e-12)
            << "u* != u^k at cell " << c << " (nu=0, zero pressure)";
        EXPECT_NEAR(uStar[c].y(), 0.0, 1e-12)
            << "v* != 0 at cell " << c << " (nu=0, zero pressure)";
    }
}

// ── MomentumSolver: pressure gradient drives momentum ───────────────────────

TEST(MomentumSolver, Solve_PressureGradient_ModifiesVelocity)
{
    // A non-uniform pressure field produces a non-zero gradient → u* ≠ u^k.
    // Just checks that the solve returns a modified field; not checking exact values.
    Mesh mesh = make4x4Mesh();
    MomentumSolver ms(mesh);

    const auto             uOld = zeroVelocity(mesh);
    Field<Eigen::Vector2d> uStar = zeroVelocity(mesh);

    // Linear pressure increasing in x: p(i,j) = x_cell_centre
    Field<double> pressure(mesh, 0.0);
    const int Nx = mesh.Nx();
    const int Ny = mesh.Ny();
    for (int i = 0; i < Nx; ++i)
        for (int j = 0; j < Ny; ++j)
            pressure(i, j) = mesh.getCellCenter(i, j).x();

    ms.solve(uStar, uOld, pressure, 0.01, 1.0, 0.01, ConvectionScheme::UPWIND);

    // With dp/dx > 0 and (1/ρ)∇p in the source, interior cells should get
    // a negative ux perturbation (pressure pushes fluid in -x direction).
    bool anyChanged = false;
    const int N = mesh.numCells();
    for (int c = 0; c < N; ++c)
        if (std::abs(uStar[c].x()) > 1e-14) { anyChanged = true; break; }

    EXPECT_TRUE(anyChanged)
        << "u* unchanged despite non-zero pressure gradient";
}

// ── MomentumSolver: invalid parameter rejection ───────────────────────────────

TEST(MomentumSolver, Solve_DtZero_Throws)
{
    Mesh mesh = make4x4Mesh();
    MomentumSolver ms(mesh);
    Field<Eigen::Vector2d> uStar = zeroVelocity(mesh);
    const auto uOld              = zeroVelocity(mesh);
    const auto pressure          = zeroPressure(mesh);

    EXPECT_THROW(ms.solve(uStar, uOld, pressure, 0.0, 1.0, 0.01, ConvectionScheme::UPWIND),
                 std::invalid_argument);
}

TEST(MomentumSolver, Solve_AlphaUZero_Throws)
{
    Mesh mesh = make4x4Mesh();
    MomentumSolver ms(mesh);
    Field<Eigen::Vector2d> uStar = zeroVelocity(mesh);
    const auto uOld              = zeroVelocity(mesh);
    const auto pressure          = zeroPressure(mesh);

    EXPECT_THROW(ms.solve(uStar, uOld, pressure, 0.01, 1.0, 0.01,
                          ConvectionScheme::UPWIND, /*alphaU=*/0.0),
                 std::invalid_argument);
}

TEST(MomentumSolver, Solve_AlphaURelaxed_ProducesFiniteResult)
{
    // With alphaU < 1 the under-relaxation modifies the diagonal (a_P → a_P/αu).
    // The two solves must give finite but different results.
    //
    // Degenerate case to avoid: uniform u^k with zero pressure gives u*=u^k
    // regardless of αu, because the under-relaxation adds zero (u^k terms cancel).
    // Instead: zero u^k with a non-zero pressure gradient forces a non-trivial
    // momentum source, so the relaxation factor changes the magnitude of u*.
    Mesh mesh = make4x4Mesh();
    MomentumSolver ms(mesh);

    const auto uOld = zeroVelocity(mesh);  // u^k = 0 → relaxation adds zero to RHS

    // Linear pressure increasing in x: produces a non-zero ∇p source.
    Field<double> pressure(mesh, 0.0);
    const int Nx = mesh.Nx();
    const int Ny = mesh.Ny();
    for (int i = 0; i < Nx; ++i)
        for (int j = 0; j < Ny; ++j)
            pressure(i, j) = mesh.getCellCenter(i, j).x();

    Field<Eigen::Vector2d> uStarFull(mesh, Eigen::Vector2d::Zero());
    Field<Eigen::Vector2d> uStarRelaxed(mesh, Eigen::Vector2d::Zero());

    ms.solve(uStarFull,    uOld, pressure, 0.01, 1.0, 0.01, ConvectionScheme::UPWIND, 1.0);
    ms.solve(uStarRelaxed, uOld, pressure, 0.01, 1.0, 0.01, ConvectionScheme::UPWIND, 0.5);

    // With αu=0.5 the diagonal doubles; u* is damped relative to αu=1 when the
    // RHS (pressure source) is unchanged (u^k=0 so no extra relaxation term).
    const int N = mesh.numCells();
    bool anyDiff = false;
    for (int c = 0; c < N; ++c)
    {
        EXPECT_TRUE(std::isfinite(uStarRelaxed[c].x()));
        EXPECT_TRUE(std::isfinite(uStarRelaxed[c].y()));
        if (std::abs(uStarRelaxed[c].x() - uStarFull[c].x()) > 1e-12)
            anyDiff = true;
    }
    EXPECT_TRUE(anyDiff) << "Relaxed and unrelaxed solves are identical — relaxation has no effect";
}

TEST(MomentumSolver, Solve_RhoZero_Throws)
{
    Mesh mesh = make4x4Mesh();
    MomentumSolver ms(mesh);
    Field<Eigen::Vector2d> uStar = zeroVelocity(mesh);
    const auto uOld              = zeroVelocity(mesh);
    const auto pressure          = zeroPressure(mesh);

    EXPECT_THROW(ms.solve(uStar, uOld, pressure, 0.01, 0.0, 0.01, ConvectionScheme::UPWIND),
                 std::invalid_argument);
}

// ── Discretization::divergenceRhieChow: sanity tests ─────────────────────────

TEST(DivergenceRhieChow, ZeroVelocityZeroPressure_ReturnsZero)
{
    // With u=0 and p=0: ū_f = 0, compact grad = 0, RC correction = 0.
    // Result must be identically zero.
    Mesh mesh = make4x4Mesh();
    const auto vel      = zeroVelocity(mesh);
    const auto pressure = zeroPressure(mesh);

    const Field<double> div =
        Discretization::divergenceRhieChow(vel, pressure, mesh, 0.01, 1.0);

    const int N = mesh.numCells();
    for (int c = 0; c < N; ++c)
        EXPECT_NEAR(div[c], 0.0, 1e-14) << "Non-zero RC divergence at cell " << c;
}

TEST(DivergenceRhieChow, UniformPressure_RCCorrectionIsZero_MatchesStandardDivergence)
{
    // For uniform pressure (p = constant), the compact face gradient and the
    // face-averaged cell-centre gradient are both zero, so the RC correction is
    // identically zero for every interior face.  divergenceRhieChow must equal
    // plain divergence regardless of the velocity field.
    //
    // Note: a *linear* pressure field does NOT give zero RC correction when the
    // gradient is computed with zero-gradient closure at boundary cells (the
    // cell-centre gradient at boundary-adjacent cells is O(h) different from the
    // exact value), which is expected and correct behaviour of the RC scheme.
    Mesh mesh = make4x4Mesh();

    // Non-zero, divergence-free velocity: u = (1, 0) everywhere.
    Field<Eigen::Vector2d> vel(mesh, Eigen::Vector2d(1.0, 0.0));

    // Uniform pressure (non-zero constant to exercise the formula non-trivially).
    const Field<double> pressure(mesh, 5.0);

    const Field<double> divRC  = Discretization::divergenceRhieChow(vel, pressure, mesh, 0.01, 1.0);
    const Field<double> divStd = Discretization::divergence(vel, mesh);

    const int N = mesh.numCells();
    for (int c = 0; c < N; ++c)
        EXPECT_NEAR(divRC[c], divStd[c], 1e-12)
            << "RC and standard divergence differ for uniform pressure at cell " << c;
}

TEST(DivergenceRhieChow, CheckerboardPressure_DiffersFromStandardDivergence)
{
    // For a checkerboard pressure field, the compact face gradient sees the
    // oscillation but the interpolated cell-centre gradient does not.
    // The RC correction is therefore non-zero, and the two divergences differ.
    Mesh mesh = make4x4Mesh();
    const auto vel = zeroVelocity(mesh);

    // Checkerboard: p(i,j) = (-1)^(i+j)
    Field<double> pressure(mesh, 0.0);
    const int Nx = mesh.Nx();
    const int Ny = mesh.Ny();
    for (int i = 0; i < Nx; ++i)
        for (int j = 0; j < Ny; ++j)
            pressure(i, j) = ((i + j) % 2 == 0) ? 1.0 : -1.0;

    const Field<double> divRC  = Discretization::divergenceRhieChow(vel, pressure, mesh, 0.01, 1.0);
    const Field<double> divStd = Discretization::divergence(vel, mesh);  // all zero (u=0)

    // Standard divergence is zero everywhere (velocity is zero).
    // RC divergence is non-zero because the checkerboard pressure has a non-zero
    // compact-minus-interpolated correction.
    bool rcNonZero = false;
    const int N = mesh.numCells();
    for (int c = 0; c < N; ++c)
    {
        EXPECT_NEAR(divStd[c], 0.0, 1e-14);
        if (std::abs(divRC[c]) > 1e-10) rcNonZero = true;
    }

    EXPECT_TRUE(rcNonZero)
        << "RC divergence is zero for checkerboard pressure — correction is missing";
}

TEST(DivergenceRhieChow, FiniteResult_NoNaN)
{
    // Arbitrary non-trivial velocity and pressure: result must be finite.
    Mesh mesh = make4x4Mesh();
    const int Nx = mesh.Nx();
    const int Ny = mesh.Ny();

    Field<Eigen::Vector2d> vel(mesh, Eigen::Vector2d::Zero());
    Field<double>          pressure(mesh, 0.0);

    // Non-uniform velocity and pressure
    for (int i = 0; i < Nx; ++i)
        for (int j = 0; j < Ny; ++j)
        {
            const double x = mesh.getCellCenter(i, j).x();
            const double y = mesh.getCellCenter(i, j).y();
            vel(i, j)      = Eigen::Vector2d(std::sin(x), std::cos(y));
            pressure(i, j) = x * x - y * y;
        }

    // Capture result via lambda to use ASSERT_NO_THROW without a default-constructed Field.
    Field<double> div(mesh, 0.0);
    ASSERT_NO_THROW(
        div = Discretization::divergenceRhieChow(vel, pressure, mesh, 0.001, 1.0));

    const int N = mesh.numCells();
    for (int c = 0; c < N; ++c)
        EXPECT_TRUE(std::isfinite(div[c])) << "Non-finite RC divergence at cell " << c;
}

// ── Boundary faces (setBoundaryFaces) ─────────────────────────────────────────

namespace
{
/// 4×8 mesh on [0, 1] × [0, 1] with a stationary bottom wall and a top wall
/// moving at U = 1; left/right unpatched (zero-gradient).
struct CouetteSetup
{
    Mesh              mesh;
    BoundaryCondition bc;

    CouetteSetup()
    {
        mesh.load(4, 8, 1.0, 1.0);
        bc.addPatch("bottom", {BoundaryType::WALL, Eigen::Vector2d(0.0, 0.0)});
        bc.addPatch("top",    {BoundaryType::WALL, Eigen::Vector2d(1.0, 0.0)});
    }

    /// Exact Couette profile u_x = y at the cell centres.
    [[nodiscard]] Field<Eigen::Vector2d> exactProfile() const
    {
        Field<Eigen::Vector2d> u(mesh, Eigen::Vector2d::Zero());
        for (int i = 0; i < mesh.Nx(); ++i)
            for (int j = 0; j < mesh.Ny(); ++j)
                u(i, j) = Eigen::Vector2d(mesh.getCellCenter(i, j).y(), 0.0);
        return u;
    }
};
}  // namespace

TEST(MomentumSolver, SetBoundaryFaces_Couette_ExactProfileIsSteady)
{
    // u_x = y satisfies the discrete viscous operator exactly when the walls
    // sit on the boundary faces (half-cell distance), so the predictor must
    // return it unchanged. The legacy zero-gradient closure would not.
    CouetteSetup s;
    MomentumSolver solver(s.mesh);
    solver.setBoundaryFaces(s.bc.buildVelocityFaces(s.mesh), s.bc.buildPressureFaces(s.mesh));

    const Field<Eigen::Vector2d> uOld = s.exactProfile();
    Field<Eigen::Vector2d>       uStar(s.mesh, Eigen::Vector2d::Zero());
    const Field<double>          p(s.mesh, 0.0);

    solver.solve(uStar, uOld, p, 1.0, 1.0, 0.1, ConvectionScheme::UPWIND, 1.0);

    for (int c = 0; c < s.mesh.numCells(); ++c)
    {
        EXPECT_NEAR(uStar[c].x(), uOld[c].x(), 1e-12) << "cell " << c;
        EXPECT_NEAR(uStar[c].y(), 0.0,         1e-12) << "cell " << c;
    }
}

TEST(MomentumSolver, SetBoundaryFaces_CouetteVariableNu_ExactProfileIsSteady)
{
    CouetteSetup s;
    MomentumSolver solver(s.mesh);
    solver.setBoundaryFaces(s.bc.buildVelocityFaces(s.mesh), s.bc.buildPressureFaces(s.mesh));

    const Field<Eigen::Vector2d> uOld = s.exactProfile();
    Field<Eigen::Vector2d>       uStar(s.mesh, Eigen::Vector2d::Zero());
    const Field<double>          p(s.mesh, 0.0);
    const Field<double>          nuEff(s.mesh, 0.1);

    solver.solve(uStar, uOld, p, nuEff, 1.0, 1.0, ConvectionScheme::CENTRAL, 1.0);

    for (int c = 0; c < s.mesh.numCells(); ++c)
        EXPECT_NEAR(uStar[c].x(), uOld[c].x(), 1e-12) << "cell " << c;
}

TEST(MomentumSolver, SetBoundaryFaces_StillFluidMovingLid_DragsTopRow)
{
    // Fluid at rest under a moving lid: only the top row feels the wall shear
    // in one implicit step's source, and it must move in +x.
    CouetteSetup s;
    MomentumSolver solver(s.mesh);
    solver.setBoundaryFaces(s.bc.buildVelocityFaces(s.mesh), s.bc.buildPressureFaces(s.mesh));

    const Field<Eigen::Vector2d> uOld(s.mesh, Eigen::Vector2d::Zero());
    Field<Eigen::Vector2d>       uStar(s.mesh, Eigen::Vector2d::Zero());
    const Field<double>          p(s.mesh, 0.0);

    solver.solve(uStar, uOld, p, 0.01, 1.0, 0.1, ConvectionScheme::UPWIND, 1.0);

    const int top = s.mesh.Ny() - 1;
    for (int i = 0; i < s.mesh.Nx(); ++i)
    {
        EXPECT_GT(uStar(i, top).x(), 0.0)              << "top row i=" << i;
        EXPECT_GT(uStar(i, top).x(), uStar(i, 0).x())  << "column i=" << i;
    }
}

TEST(MomentumSolver, SetBoundaryFaces_MismatchedLists_Throws)
{
    CouetteSetup s;
    MomentumSolver solver(s.mesh);
    std::vector<PressureBoundaryFace> pf = s.bc.buildPressureFaces(s.mesh);
    pf.pop_back();

    EXPECT_THROW(solver.setBoundaryFaces(s.bc.buildVelocityFaces(s.mesh), std::move(pf)),
                 std::invalid_argument);
}
