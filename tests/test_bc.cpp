#include "core/BoundaryCondition.hpp"
#include "core/Field.hpp"
#include "core/Mesh.hpp"

#include <Eigen/Dense>
#include <gtest/gtest.h>
#include <stdexcept>

namespace
{
Mesh makeSquareMesh()
{
    Mesh m;
    m.load(8, 8, 1.0, 1.0);
    return m;
}
}  // namespace

// ── addPatch / getPatch ───────────────────────────────────────────────────────

TEST(BoundaryCondition, AddAndGetPatch_RoundTrip)
{
    BoundaryCondition bc;
    BoundaryPatch patch{BoundaryType::WALL, Eigen::Vector2d(1.0, 0.0)};
    bc.addPatch("top", patch);

    const BoundaryPatch& retrieved = bc.getPatch("top");
    EXPECT_EQ(retrieved.type, BoundaryType::WALL);
    EXPECT_NEAR(retrieved.value.x(), 1.0, 1e-12);
    EXPECT_NEAR(retrieved.value.y(), 0.0, 1e-12);
}

TEST(BoundaryCondition, GetPatch_UnknownName_Throws)
{
    BoundaryCondition bc;
    EXPECT_THROW(static_cast<void>(bc.getPatch("nonexistent")), std::out_of_range);
}

// ── applyVelocity ─────────────────────────────────────────────────────────────

/// Lid-driven cavity setup: top WALL at (1,0), all other walls at (0,0).
/// After applyVelocity:
///   - top row (j=7): all cells = (1,0)  — top applied last, wins at corners
///   - bottom row (j=0): all cells = (0,0)
///   - left/right columns (j<7): = (0,0)
TEST(BoundaryCondition, ApplyVelocity_LidWall_SetsLidVelocity)
{
    Mesh mesh = makeSquareMesh();
    Field<Eigen::Vector2d> vel(mesh, Eigen::Vector2d::Zero());

    BoundaryCondition bc;
    bc.addPatch("top",    {BoundaryType::WALL, Eigen::Vector2d(1.0, 0.0)});
    bc.addPatch("bottom", {BoundaryType::WALL, Eigen::Vector2d::Zero()});
    bc.addPatch("left",   {BoundaryType::WALL, Eigen::Vector2d::Zero()});
    bc.addPatch("right",  {BoundaryType::WALL, Eigen::Vector2d::Zero()});

    EXPECT_NO_THROW(bc.applyVelocity(vel, mesh));

    const int Nx = 8, Ny = 8;

    // Top row: all i at j=Ny-1 must equal lid velocity (1,0).
    // Top is applied last, so it wins at the top corners even if left/right
    // wrote (0,0) there first.
    for (int i = 0; i < Nx; ++i) {
        EXPECT_NEAR(vel(i, Ny - 1).x(), 1.0, 1e-12) << "top row, i=" << i;
        EXPECT_NEAR(vel(i, Ny - 1).y(), 0.0, 1e-12) << "top row, i=" << i;
    }

    // Bottom row: all (0,0).
    for (int i = 0; i < Nx; ++i) {
        EXPECT_NEAR(vel(i, 0).x(), 0.0, 1e-12) << "bottom row, i=" << i;
        EXPECT_NEAR(vel(i, 0).y(), 0.0, 1e-12) << "bottom row, i=" << i;
    }

    // Left column, rows below the top corner: (0,0).
    for (int j = 0; j < Ny - 1; ++j) {
        EXPECT_NEAR(vel(0, j).x(), 0.0, 1e-12) << "left col, j=" << j;
        EXPECT_NEAR(vel(0, j).y(), 0.0, 1e-12) << "left col, j=" << j;
    }

    // Right column, rows below the top corner: (0,0).
    for (int j = 0; j < Ny - 1; ++j) {
        EXPECT_NEAR(vel(Nx - 1, j).x(), 0.0, 1e-12) << "right col, j=" << j;
        EXPECT_NEAR(vel(Nx - 1, j).y(), 0.0, 1e-12) << "right col, j=" << j;
    }
}

/// OUTLET patch: zero-gradient policy → boundary cells must be left unchanged.
TEST(BoundaryCondition, ApplyVelocity_Outlet_LeavesFieldUnchanged)
{
    Mesh mesh = makeSquareMesh();
    const Eigen::Vector2d initVel(2.0, 3.0);
    Field<Eigen::Vector2d> vel(mesh, initVel);

    BoundaryCondition bc;
    bc.addPatch("right", {BoundaryType::OUTLET, Eigen::Vector2d::Zero()});

    bc.applyVelocity(vel, mesh);

    const int Nx = 8, Ny = 8;
    // Right column must be unchanged (zero-gradient = leave as-is).
    for (int j = 0; j < Ny; ++j) {
        EXPECT_NEAR(vel(Nx - 1, j).x(), initVel.x(), 1e-12) << "right col, j=" << j;
        EXPECT_NEAR(vel(Nx - 1, j).y(), initVel.y(), 1e-12) << "right col, j=" << j;
    }
}

/// Field from different mesh must throw std::invalid_argument.
TEST(BoundaryCondition, ApplyVelocity_MeshMismatch_Throws)
{
    Mesh mesh1 = makeSquareMesh();   // 8x8
    Mesh mesh2;
    mesh2.load(4, 4, 1.0, 1.0);     // 4x4

    Field<Eigen::Vector2d> vel(mesh1, Eigen::Vector2d::Zero());

    BoundaryCondition bc;
    bc.addPatch("top", {BoundaryType::WALL, Eigen::Vector2d(1.0, 0.0)});

    EXPECT_THROW(bc.applyVelocity(vel, mesh2), std::invalid_argument);
}

// ── applyPressure ─────────────────────────────────────────────────────────────

/// OUTLET patch: right boundary cells must be set to 0; interior unchanged.
TEST(BoundaryCondition, ApplyPressure_Outlet_ZeroPressureBoundary)
{
    Mesh mesh = makeSquareMesh();
    Field<double> p(mesh, 1.0);  // start with non-zero pressure everywhere

    BoundaryCondition bc;
    bc.addPatch("right", {BoundaryType::OUTLET, Eigen::Vector2d::Zero()});

    EXPECT_NO_THROW(bc.applyPressure(p, mesh));

    const int Nx = 8, Ny = 8;

    // Right boundary column: must all be 0.
    for (int j = 0; j < Ny; ++j) {
        EXPECT_NEAR(p(Nx - 1, j), 0.0, 1e-12) << "right col, j=" << j;
    }

    // An interior cell must remain unchanged at 1.0.
    EXPECT_NEAR(p(Nx / 2, Ny / 2), 1.0, 1e-12);
}

/// WALL patch on pressure: zero-gradient → boundary cells must be left unchanged.
TEST(BoundaryCondition, ApplyPressure_Wall_LeavesFieldUnchanged)
{
    Mesh mesh = makeSquareMesh();
    const double initP = 5.0;
    Field<double> p(mesh, initP);

    BoundaryCondition bc;
    bc.addPatch("top", {BoundaryType::WALL, Eigen::Vector2d::Zero()});

    bc.applyPressure(p, mesh);

    const int Nx = 8, Ny = 8;
    // Top row must be unchanged.
    for (int i = 0; i < Nx; ++i) {
        EXPECT_NEAR(p(i, Ny - 1), initP, 1e-12) << "top row, i=" << i;
    }
}

/// Field from different mesh must throw std::invalid_argument.
TEST(BoundaryCondition, ApplyPressure_MeshMismatch_Throws)
{
    Mesh mesh1 = makeSquareMesh();   // 8x8
    Mesh mesh2;
    mesh2.load(4, 4, 1.0, 1.0);     // 4x4

    Field<double> p(mesh1, 0.0);

    BoundaryCondition bc;
    bc.addPatch("right", {BoundaryType::OUTLET, Eigen::Vector2d::Zero()});

    EXPECT_THROW(bc.applyPressure(p, mesh2), std::invalid_argument);
}
