#include "mesh/UnstructuredField.hpp"
#include "mesh/UnstructuredMesh.hpp"

#include <gtest/gtest.h>

#include <array>
#include <cmath>

// ── Fixture: simple 2-cell mesh ──────────────────────────────────────────────

/// Builds a minimal 2-quad mesh (2×1 grid) for field tests.
/// 6 nodes, 2 quad cells, no boundary patches registered.
static UnstructuredMesh buildSimpleMesh()
{
    UnstructuredMesh mesh;

    //  3---4---5
    //  |   |   |
    //  0---1---2
    mesh.addNode({{0.0, 0.0}});  // 0
    mesh.addNode({{1.0, 0.0}});  // 1
    mesh.addNode({{2.0, 0.0}});  // 2
    mesh.addNode({{0.0, 1.0}});  // 3
    mesh.addNode({{1.0, 1.0}});  // 4
    mesh.addNode({{2.0, 1.0}});  // 5

    // Left quad: CCW 0,1,4,3
    std::array<NodeId, 4> quadL{NodeId{0}, NodeId{1}, NodeId{4}, NodeId{3}};
    // Right quad: CCW 1,2,5,4
    std::array<NodeId, 4> quadR{NodeId{1}, NodeId{2}, NodeId{5}, NodeId{4}};

    mesh.addCell(std::span<const NodeId>{quadL.data(), 4});
    mesh.addCell(std::span<const NodeId>{quadR.data(), 4});

    // Build connectivity with no explicit boundary patches
    // (all edges become boundary faces with patchId=-1).
    std::unordered_map<std::pair<int, int>, int, PairHash> bdry;
    mesh.buildConnectivity(bdry);
    return mesh;
}

// ── UnstructuredField<double> tests ──────────────────────────────────────────

TEST(UnstructuredField, Construction_SizeMatchesMesh)
{
    const UnstructuredMesh mesh = buildSimpleMesh();
    const UnstructuredField<double> f(mesh);
    EXPECT_EQ(f.size(), mesh.numCells());
}

TEST(UnstructuredField, Construction_WithInitialValue)
{
    const UnstructuredMesh mesh = buildSimpleMesh();
    const UnstructuredField<double> f(mesh, 3.14);
    for (int ci = 0; ci < mesh.numCells(); ++ci)
        EXPECT_NEAR(f[CellId{ci}], 3.14, 1e-15);
}

TEST(UnstructuredField, SetAll)
{
    const UnstructuredMesh mesh = buildSimpleMesh();
    UnstructuredField<double> f(mesh, 0.0);
    f.setAll(7.0);
    for (int ci = 0; ci < mesh.numCells(); ++ci)
        EXPECT_NEAR(f[CellId{ci}], 7.0, 1e-15);
}

TEST(UnstructuredField, ReadWriteByIndex)
{
    const UnstructuredMesh mesh = buildSimpleMesh();
    UnstructuredField<double> f(mesh, 0.0);
    f[CellId{0}] = 1.0;
    f[CellId{1}] = 2.0;
    EXPECT_NEAR(f[CellId{0}], 1.0, 1e-15);
    EXPECT_NEAR(f[CellId{1}], 2.0, 1e-15);
}

TEST(UnstructuredField, Norm_AllOnes)
{
    // norm of N all-1 values = sqrt(N)
    const UnstructuredMesh mesh = buildSimpleMesh();
    const UnstructuredField<double> f(mesh, 1.0);
    const double expected = std::sqrt(static_cast<double>(mesh.numCells()));
    EXPECT_NEAR(f.norm(), expected, 1e-12);
}

TEST(UnstructuredField, Norm_AllZero)
{
    const UnstructuredMesh mesh = buildSimpleMesh();
    const UnstructuredField<double> f(mesh, 0.0);
    EXPECT_NEAR(f.norm(), 0.0, 1e-15);
}

TEST(UnstructuredField, OperatorPlus)
{
    const UnstructuredMesh mesh = buildSimpleMesh();
    UnstructuredField<double> a(mesh, 2.0);
    UnstructuredField<double> b(mesh, 3.0);
    const UnstructuredField<double> c = a + b;
    for (int ci = 0; ci < mesh.numCells(); ++ci)
        EXPECT_NEAR(c[CellId{ci}], 5.0, 1e-15);
}

TEST(UnstructuredField, OperatorScalar)
{
    const UnstructuredMesh mesh = buildSimpleMesh();
    const UnstructuredField<double> a(mesh, 4.0);
    const UnstructuredField<double> b = a * 2.5;
    for (int ci = 0; ci < mesh.numCells(); ++ci)
        EXPECT_NEAR(b[CellId{ci}], 10.0, 1e-15);
}

TEST(UnstructuredField, MeshReference)
{
    const UnstructuredMesh mesh = buildSimpleMesh();
    const UnstructuredField<double> f(mesh);
    EXPECT_EQ(&f.mesh(), &mesh);
}

// ── UnstructuredField<Eigen::Vector2d> tests ─────────────────────────────────

TEST(UnstructuredFieldVector, Construction_WithInitialValue)
{
    const UnstructuredMesh mesh = buildSimpleMesh();
    const Eigen::Vector2d v{1.0, 2.0};
    const UnstructuredField<Eigen::Vector2d> f(mesh, v);
    for (int ci = 0; ci < mesh.numCells(); ++ci)
    {
        EXPECT_NEAR(f[CellId{ci}].x(), 1.0, 1e-15);
        EXPECT_NEAR(f[CellId{ci}].y(), 2.0, 1e-15);
    }
}

TEST(UnstructuredFieldVector, Norm_AllUnitVectors)
{
    // 2 cells, each with vector (1,0) → norm = sqrt(2*1²) = sqrt(2)
    const UnstructuredMesh mesh = buildSimpleMesh();
    const UnstructuredField<Eigen::Vector2d> f(mesh, Eigen::Vector2d{1.0, 0.0});
    const double expected = std::sqrt(static_cast<double>(mesh.numCells()));
    EXPECT_NEAR(f.norm(), expected, 1e-12);
}

TEST(UnstructuredFieldVector, OperatorPlus)
{
    const UnstructuredMesh mesh = buildSimpleMesh();
    const UnstructuredField<Eigen::Vector2d> a(mesh, Eigen::Vector2d{1.0, 0.0});
    const UnstructuredField<Eigen::Vector2d> b(mesh, Eigen::Vector2d{0.0, 1.0});
    const UnstructuredField<Eigen::Vector2d> c = a + b;
    for (int ci = 0; ci < mesh.numCells(); ++ci)
    {
        EXPECT_NEAR(c[CellId{ci}].x(), 1.0, 1e-15);
        EXPECT_NEAR(c[CellId{ci}].y(), 1.0, 1e-15);
    }
}

TEST(UnstructuredFieldVector, OperatorScalar)
{
    const UnstructuredMesh mesh = buildSimpleMesh();
    const UnstructuredField<Eigen::Vector2d> a(mesh, Eigen::Vector2d{3.0, 4.0});
    const UnstructuredField<Eigen::Vector2d> b = a * 2.0;
    for (int ci = 0; ci < mesh.numCells(); ++ci)
    {
        EXPECT_NEAR(b[CellId{ci}].x(), 6.0, 1e-15);
        EXPECT_NEAR(b[CellId{ci}].y(), 8.0, 1e-15);
    }
}
