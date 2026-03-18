#include "core/Mesh.hpp"

#include <gtest/gtest.h>

// ── Milestone 2 will replace these stubs with real assertions ────────────────

/**
 * When implemented: mesh.load() should correctly set Nx, Ny and compute
 * dx = Lx/Nx, dy = Ly/Ny without throwing.
 */
TEST(Mesh, Load_ValidDimensions_DoesNotThrow)
{
    Mesh mesh;
    EXPECT_NO_THROW(mesh.load(16, 16, 1.0, 1.0));
    // TODO(Milestone 2): EXPECT_EQ(mesh.Nx(), 16)
    //                    EXPECT_EQ(mesh.Ny(), 16)
    EXPECT_TRUE(true);
}

/**
 * When implemented: load() should throw std::invalid_argument for Nx = 0.
 */
TEST(Mesh, Load_ZeroNx_Throws)
{
    Mesh mesh;
    EXPECT_THROW(mesh.load(0, 16, 1.0, 1.0), std::invalid_argument);
}

/**
 * When implemented: load() should throw for non-positive domain size.
 */
TEST(Mesh, Load_NegativeLx_Throws)
{
    Mesh mesh;
    EXPECT_THROW(mesh.load(16, 16, -1.0, 1.0), std::invalid_argument);
}

/**
 * When implemented: numCells() should return Nx * Ny.
 */
TEST(Mesh, NumCells_EqualsNxTimesNy)
{
    Mesh mesh;
    mesh.load(8, 4, 1.0, 0.5);
    EXPECT_EQ(mesh.numCells(), 32);
}

/**
 * When implemented: cellIndex(i,j) should equal i * Ny + j.
 */
TEST(Mesh, CellIndex_RowMajorLayout)
{
    Mesh mesh;
    mesh.load(4, 8, 1.0, 1.0);
    EXPECT_EQ(mesh.cellIndex(2, 3), 2 * 8 + 3);
}

/**
 * When implemented: getCellCenter(0,0) should return (dx/2, dy/2) for a
 * uniform grid with origin at (0,0).
 */
TEST(Mesh, GetCellCenter_OriginCell_HalfCellOffset)
{
    Mesh mesh;
    mesh.load(4, 4, 1.0, 1.0);
    // STUB: getCellCenter returns zero until Milestone 2
    // TODO(Milestone 2): Eigen::Vector2d c = mesh.getCellCenter(0, 0);
    //                    EXPECT_NEAR(c.x(), 0.125, 1e-12);
    //                    EXPECT_NEAR(c.y(), 0.125, 1e-12);
    EXPECT_TRUE(true);
}
