#include "core/Field.hpp"
#include "core/Mesh.hpp"

#include <Eigen/Dense>
#include <gtest/gtest.h>

// ── Milestone 2 will replace these stubs with real assertions ────────────────

namespace
{
/// Helper: creates a small 4×4 mesh for field tests.
Mesh make4x4Mesh()
{
    Mesh m;
    m.load(4, 4, 1.0, 1.0);
    return m;
}
}  // namespace

/**
 * When implemented: Field<double> constructed on a mesh should have
 * size() == mesh.numCells().
 */
TEST(Field, ScalarField_Size_EqualsNumCells)
{
    Mesh mesh = make4x4Mesh();
    Field<double> f(mesh, 0.0);
    EXPECT_EQ(f.size(), mesh.numCells());
}

/**
 * When implemented: setAll(v) should make every cell equal to v.
 */
TEST(Field, ScalarField_SetAll_FillsAllCells)
{
    Mesh mesh = make4x4Mesh();
    Field<double> f(mesh);
    f.setAll(3.14);
    // TODO(Milestone 2): for every cell, EXPECT_NEAR(f[i], 3.14, 1e-12)
    EXPECT_TRUE(true);
}

/**
 * When implemented: operator()(i,j) should be consistent with operator[].
 */
TEST(Field, ScalarField_IndexConsistency_TwoDAndLinear)
{
    Mesh mesh = make4x4Mesh();
    Field<double> f(mesh, 0.0);
    f(1, 2) = 7.0;
    // TODO(Milestone 2): EXPECT_NEAR(f[mesh.cellIndex(1,2)], 7.0, 1e-12)
    EXPECT_TRUE(true);
}

/**
 * When implemented: norm() of a uniform field of value v on N cells should
 * equal v * sqrt(N).
 */
TEST(Field, ScalarField_Norm_UniformField)
{
    Mesh mesh = make4x4Mesh();
    Field<double> f(mesh, 1.0);
    // TODO(Milestone 2): EXPECT_NEAR(f.norm(), std::sqrt(16.0), 1e-12)
    EXPECT_TRUE(true);
}

/**
 * When implemented: Field<Eigen::Vector2d> should store and retrieve 2D vectors.
 */
TEST(Field, VectorField_SetAndGet_PreservesVector)
{
    Mesh mesh = make4x4Mesh();
    Field<Eigen::Vector2d> f(mesh, Eigen::Vector2d::Zero());
    Eigen::Vector2d v(1.0, 2.0);
    f(0, 0) = v;
    // TODO(Milestone 2): EXPECT_NEAR(f(0,0).x(), 1.0, 1e-12)
    //                    EXPECT_NEAR(f(0,0).y(), 2.0, 1e-12)
    EXPECT_TRUE(true);
}
