#include "solver/WallDistance.hpp"
#include "core/BoundaryCondition.hpp"
#include "core/Mesh.hpp"

#include <gtest/gtest.h>
#include <cmath>

// ── Helpers ───────────────────────────────────────────────────────────────────

namespace
{

/// Builds a 4×4 unit-square mesh (Lx=1, Ly=1).
Mesh make4x4Mesh()
{
    Mesh m;
    m.load(4, 4, 1.0, 1.0);
    return m;
}

/// Returns a BoundaryCondition with only the bottom edge as WALL.
BoundaryCondition bottomWallBC()
{
    BoundaryCondition bc;
    bc.addPatch("bottom", {BoundaryType::WALL, Eigen::Vector2d::Zero()});
    return bc;
}

} // namespace

// ── WallDistance: monotonic increase away from bottom wall ───────────────────

TEST(WallDistance, BottomWall_MonotonicIncrease)
{
    // On a 4×4 mesh with bottom=WALL (j=0 row), wall distance should increase
    // monotonically as j increases.  All cells in the same column (same i) should
    // have the same distance because the mesh is uniform.
    const Mesh             mesh = make4x4Mesh();
    const BoundaryCondition bc  = bottomWallBC();

    const Field<double> dist = WallDistance::compute(mesh, bc);

    // dy = 1.0 / 4 = 0.25
    // j=0: dist = 0.5 * dy = 0.125
    // j=1: dist = 1.5 * dy = 0.375
    // j=2: dist = 2.5 * dy = 0.625
    // j=3: dist = 3.5 * dy = 0.875
    const double dy = 1.0 / 4.0;

    for (int i = 0; i < mesh.Nx(); ++i)
    {
        for (int j = 0; j < mesh.Ny(); ++j)
        {
            const double expected = (static_cast<double>(j) + 0.5) * dy;
            EXPECT_NEAR(dist(i, j), expected, 1.0e-10)
                << "Incorrect wall distance at (" << i << ", " << j << ")";
        }
    }
}

TEST(WallDistance, BottomWall_MonotonicByRow)
{
    // Wall distances must increase monotonically as j increases.
    const Mesh             mesh = make4x4Mesh();
    const BoundaryCondition bc  = bottomWallBC();

    const Field<double> dist = WallDistance::compute(mesh, bc);

    for (int i = 0; i < mesh.Nx(); ++i)
    {
        for (int j = 1; j < mesh.Ny(); ++j)
        {
            EXPECT_GT(dist(i, j), dist(i, j - 1))
                << "Wall distance not increasing at column " << i
                << " from j=" << (j-1) << " to j=" << j;
        }
    }
}

TEST(WallDistance, BottomWall_WallRowHasMinDist)
{
    // The j=0 row (wall cells) must have the smallest distance in their column.
    const Mesh             mesh = make4x4Mesh();
    const BoundaryCondition bc  = bottomWallBC();

    const Field<double> dist = WallDistance::compute(mesh, bc);

    for (int i = 0; i < mesh.Nx(); ++i)
        for (int j = 1; j < mesh.Ny(); ++j)
            EXPECT_LT(dist(i, 0), dist(i, j))
                << "Wall row j=0 does not have minimum distance in column " << i;
}

TEST(WallDistance, NoWall_DistancesRemainMax)
{
    // A BoundaryCondition with no WALL patches should leave all distances at max.
    const Mesh              mesh = make4x4Mesh();
    const BoundaryCondition bc;  // empty: no patches

    const Field<double> dist = WallDistance::compute(mesh, bc);

    const int N = mesh.numCells();
    for (int c = 0; c < N; ++c)
        EXPECT_GT(dist[c], 1.0e10)
            << "Expected max distance (no wall), got " << dist[c] << " at cell " << c;
}
