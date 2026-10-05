#include "mesh/UnstructuredMesh.hpp"

#include <gtest/gtest.h>

#include <array>
#include <cmath>

// ── Helpers ───────────────────────────────────────────────────────────────────

/// Builds a 2-triangle mesh from a unit square split along the diagonal.
/// Nodes (CCW):
///   0:(0,0), 1:(1,0), 2:(1,1), 3:(0,1)
/// Triangle A: nodes 0,1,2  (lower-right)
/// Triangle B: nodes 0,2,3  (upper-left)
/// Shared interior face: edge 0-2
/// Boundary faces: 0-1, 1-2, 2-3, 3-0  (4 boundary)
/// Total faces: 5 (1 interior + 4 boundary)
static UnstructuredMesh buildTwoTriangleMesh()
{
    UnstructuredMesh mesh;

    const NodeId n0 = mesh.addNode({{0.0, 0.0}});
    const NodeId n1 = mesh.addNode({{1.0, 0.0}});
    const NodeId n2 = mesh.addNode({{1.0, 1.0}});
    const NodeId n3 = mesh.addNode({{0.0, 1.0}});

    // Register 4 boundary patches (one per boundary edge)
    const int pBot  = mesh.addPatch("bottom");  // 0-1
    const int pRight= mesh.addPatch("right");   // 1-2
    const int pTop  = mesh.addPatch("top");     // 2-3
    const int pLeft = mesh.addPatch("left");    // 3-0

    std::array<NodeId, 3> triA{n0, n1, n2};
    std::array<NodeId, 3> triB{n0, n2, n3};

    mesh.addCell(std::span<const NodeId>{triA.data(), 3});
    mesh.addCell(std::span<const NodeId>{triB.data(), 3});

    // Build boundary edge map for all 4 boundary edges.
    std::unordered_map<std::pair<int, int>, int, PairHash> bdry;
    // bottom: node 0-1
    bdry[{0, 1}] = pBot;
    // right: node 1-2
    bdry[{1, 2}] = pRight;
    // top: node 2-3
    bdry[{2, 3}] = pTop;
    // left: node 0-3
    bdry[{0, 3}] = pLeft;

    mesh.buildConnectivity(bdry);
    return mesh;
}

// ── Tests ─────────────────────────────────────────────────────────────────────

TEST(UnstructuredMesh, NodeCount)
{
    UnstructuredMesh mesh;
    mesh.addNode({{0.0, 0.0}});
    mesh.addNode({{1.0, 0.0}});
    EXPECT_EQ(mesh.numNodes(), 2);
}

TEST(UnstructuredMesh, AddNode_ReturnsIncrementingId)
{
    UnstructuredMesh mesh;
    const NodeId id0 = mesh.addNode({{0.0, 0.0}});
    const NodeId id1 = mesh.addNode({{1.0, 0.0}});
    EXPECT_EQ(toInt(id0), 0);
    EXPECT_EQ(toInt(id1), 1);
}

TEST(UnstructuredMesh, AddPatch_IdempotentOnSameName)
{
    UnstructuredMesh mesh;
    const int idx0 = mesh.addPatch("inlet");
    const int idx1 = mesh.addPatch("inlet");
    EXPECT_EQ(idx0, idx1);
    EXPECT_EQ(mesh.numPatches(), 1);
}

TEST(UnstructuredMesh, BuildConnectivity_CellCount)
{
    const UnstructuredMesh mesh = buildTwoTriangleMesh();
    EXPECT_EQ(mesh.numCells(), 2);
}

TEST(UnstructuredMesh, BuildConnectivity_NumInteriorFaces)
{
    // Two triangles sharing one edge → 1 interior face.
    const UnstructuredMesh mesh = buildTwoTriangleMesh();
    EXPECT_EQ(mesh.numInteriorFaces(), 1);
}

TEST(UnstructuredMesh, BuildConnectivity_TotalFaces)
{
    // 2 triangles × 3 edges − 1 shared edge = 5 unique faces.
    const UnstructuredMesh mesh = buildTwoTriangleMesh();
    EXPECT_EQ(mesh.numFaces(), 5);
}

TEST(UnstructuredMesh, BuildConnectivity_BoundaryFaces)
{
    const UnstructuredMesh mesh = buildTwoTriangleMesh();
    const int nBdry = mesh.numFaces() - mesh.numInteriorFaces();
    EXPECT_EQ(nBdry, 4);
}

TEST(UnstructuredMesh, BuildConnectivity_CellCentroid_TriangleA)
{
    // Triangle A: (0,0),(1,0),(1,1) → centroid = (2/3, 1/3)
    const UnstructuredMesh mesh = buildTwoTriangleMesh();
    const MeshCell& ca = mesh.cell(CellId{0});
    EXPECT_NEAR(ca.centre.x(), 2.0 / 3.0, 1e-10);
    EXPECT_NEAR(ca.centre.y(), 1.0 / 3.0, 1e-10);
}

TEST(UnstructuredMesh, BuildConnectivity_CellCentroid_TriangleB)
{
    // Triangle B: (0,0),(1,1),(0,1) → centroid = (1/3, 2/3)
    const UnstructuredMesh mesh = buildTwoTriangleMesh();
    const MeshCell& cb = mesh.cell(CellId{1});
    EXPECT_NEAR(cb.centre.x(), 1.0 / 3.0, 1e-10);
    EXPECT_NEAR(cb.centre.y(), 2.0 / 3.0, 1e-10);
}

TEST(UnstructuredMesh, BuildConnectivity_CellVolume)
{
    // Each triangle has area = 0.5 (half of unit square)
    const UnstructuredMesh mesh = buildTwoTriangleMesh();
    EXPECT_NEAR(mesh.cell(CellId{0}).volume, 0.5, 1e-10);
    EXPECT_NEAR(mesh.cell(CellId{1}).volume, 0.5, 1e-10);
}

TEST(UnstructuredMesh, BuildConnectivity_InteriorFace_HasNeighbour)
{
    const UnstructuredMesh mesh = buildTwoTriangleMesh();
    // Interior face is at index 0.
    const MeshFace& f = mesh.face(FaceId{0});
    EXPECT_GE(toInt(f.neighbour), 0);
}

TEST(UnstructuredMesh, BuildConnectivity_BoundaryFace_NoNeighbour)
{
    const UnstructuredMesh mesh = buildTwoTriangleMesh();
    // Boundary faces start at index numInteriorFaces.
    const int bStart = mesh.numInteriorFaces();
    EXPECT_LT(bStart, mesh.numFaces());
    const MeshFace& fb = mesh.face(FaceId{bStart});
    EXPECT_LT(toInt(fb.neighbour), 0);
}

TEST(UnstructuredMesh, BuildConnectivity_FaceNormal_UnitLength)
{
    const UnstructuredMesh mesh = buildTwoTriangleMesh();
    for (int fi = 0; fi < mesh.numFaces(); ++fi)
    {
        const double nlen = mesh.face(FaceId{fi}).normal.norm();
        EXPECT_NEAR(nlen, 1.0, 1e-10) << "Face " << fi << " normal not unit length";
    }
}

TEST(UnstructuredMesh, BuildConnectivity_PatchFaces_BottomHasOneFace)
{
    // Bottom edge (0-1) is one face.
    const UnstructuredMesh mesh = buildTwoTriangleMesh();
    const auto& faces = mesh.patchFaces("bottom");
    EXPECT_EQ(static_cast<int>(faces.size()), 1);
}

TEST(UnstructuredMesh, BuildConnectivity_AllPatches_TotalBoundaryFaces)
{
    // Sum of patch face counts must equal total boundary face count.
    const UnstructuredMesh mesh = buildTwoTriangleMesh();
    int total = 0;
    for (int pi = 0; pi < mesh.numPatches(); ++pi)
        total += static_cast<int>(mesh.patchFaces(pi).size());
    EXPECT_EQ(total, mesh.numFaces() - mesh.numInteriorFaces());
}

TEST(UnstructuredMesh, IsBoundaryFace)
{
    const UnstructuredMesh mesh = buildTwoTriangleMesh();
    // First face is interior.
    EXPECT_FALSE(mesh.isBoundaryFace(FaceId{0}));
    // Faces from numInteriorFaces onward are boundary.
    EXPECT_TRUE(mesh.isBoundaryFace(FaceId{mesh.numInteriorFaces()}));
}

TEST(UnstructuredMesh, CellFaceCount_Triangle)
{
    // Each triangle has exactly 3 bounding faces.
    const UnstructuredMesh mesh = buildTwoTriangleMesh();
    EXPECT_EQ(static_cast<int>(mesh.cell(CellId{0}).faces.size()), 3);
    EXPECT_EQ(static_cast<int>(mesh.cell(CellId{1}).faces.size()), 3);
}
