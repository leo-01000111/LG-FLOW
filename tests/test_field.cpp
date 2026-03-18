#include "core/Field.hpp"
#include "core/Mesh.hpp"

#include <Eigen/Dense>
#include <cmath>
#include <gtest/gtest.h>

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

TEST(Field, ScalarField_Size_EqualsNumCells)
{
    Mesh mesh = make4x4Mesh();
    Field<double> f(mesh, 0.0);
    EXPECT_EQ(f.size(), mesh.numCells());
}

TEST(Field, ScalarField_SetAll_FillsAllCells)
{
    Mesh mesh = make4x4Mesh();
    Field<double> f(mesh);
    f.setAll(3.14);
    for (int k = 0; k < f.size(); ++k)
        EXPECT_NEAR(f[k], 3.14, 1e-12);
}

TEST(Field, ScalarField_IndexConsistency_TwoDAndLinear)
{
    Mesh mesh = make4x4Mesh();
    Field<double> f(mesh, 0.0);
    f(1, 2) = 7.0;
    EXPECT_NEAR(f[mesh.cellIndex(1, 2)], 7.0, 1e-12);
    EXPECT_NEAR(f[mesh.cellIndex(0, 0)], 0.0, 1e-12);
}

TEST(Field, ScalarField_Norm_UniformField)
{
    Mesh mesh = make4x4Mesh();
    Field<double> f(mesh, 1.0);
    // 16 cells each with value 1.0: norm = sqrt(16 * 1^2) = 4
    EXPECT_NEAR(f.norm(), std::sqrt(16.0), 1e-12);
}

TEST(Field, VectorField_Norm_UniformField)
{
    Mesh mesh = make4x4Mesh();
    // |v|^2 = 3^2 + 4^2 = 25  →  norm = sqrt(16 * 25) = 20
    Field<Eigen::Vector2d> f(mesh, Eigen::Vector2d(3.0, 4.0));
    EXPECT_NEAR(f.norm(), std::sqrt(16.0 * 25.0), 1e-12);
}

TEST(Field, ScalarField_OperatorPlus_SumIsCorrect)
{
    Mesh mesh = make4x4Mesh();
    Field<double> a(mesh, 2.0);
    Field<double> b(mesh, 3.0);
    Field<double> c = a + b;
    for (int k = 0; k < c.size(); ++k)
        EXPECT_NEAR(c[k], 5.0, 1e-12);
}

TEST(Field, ScalarField_OperatorMultiply_ScaledIsCorrect)
{
    Mesh mesh = make4x4Mesh();
    Field<double> f(mesh, 4.0);
    Field<double> g = f * 2.5;
    for (int k = 0; k < g.size(); ++k)
        EXPECT_NEAR(g[k], 10.0, 1e-12);
}

TEST(Field, ScalarField_OperatorPlus_MeshMismatch_Throws)
{
    Mesh mesh1 = make4x4Mesh();
    Mesh mesh2 = make4x4Mesh();
    Field<double> a(mesh1, 1.0);
    Field<double> b(mesh2, 1.0);
    EXPECT_THROW(static_cast<void>(a + b), std::invalid_argument);
}

TEST(Field, VectorField_SetAndGet_PreservesVector)
{
    Mesh mesh = make4x4Mesh();
    Field<Eigen::Vector2d> f(mesh, Eigen::Vector2d::Zero());
    Eigen::Vector2d v(1.0, 2.0);
    f(0, 0) = v;
    EXPECT_NEAR(f(0, 0).x(), 1.0, 1e-12);
    EXPECT_NEAR(f(0, 0).y(), 2.0, 1e-12);
    // Other cells remain zero
    EXPECT_NEAR(f(1, 0).norm(), 0.0, 1e-12);
}
