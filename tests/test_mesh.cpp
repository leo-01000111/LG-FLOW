#include "core/Mesh.hpp"

#include <gtest/gtest.h>

// ── Helpers ───────────────────────────────────────────────────────────────────

namespace
{
/// Returns a 4×4 uniform mesh over [0,1]×[0,1].
/// dx = dy = 0.25, numCells = 16, numFaces = 40.
Mesh make4x4()
{
    Mesh m;
    m.load(4, 4, 1.0, 1.0);
    return m;
}
}  // namespace

// ── Load validation ───────────────────────────────────────────────────────────

TEST(Mesh, Load_ValidDimensions_DoesNotThrow)
{
    Mesh mesh;
    EXPECT_NO_THROW(mesh.load(16, 16, 1.0, 1.0));
    EXPECT_EQ(mesh.Nx(), 16);
    EXPECT_EQ(mesh.Ny(), 16);
}

TEST(Mesh, Load_ZeroNx_Throws)
{
    Mesh mesh;
    EXPECT_THROW(mesh.load(0, 16, 1.0, 1.0), std::invalid_argument);
}

TEST(Mesh, Load_NegativeLx_Throws)
{
    Mesh mesh;
    EXPECT_THROW(mesh.load(16, 16, -1.0, 1.0), std::invalid_argument);
}

// ── Cell count and indexing ───────────────────────────────────────────────────

TEST(Mesh, NumCells_EqualsNxTimesNy)
{
    Mesh mesh;
    mesh.load(8, 4, 1.0, 0.5);
    EXPECT_EQ(mesh.numCells(), 32);
}

TEST(Mesh, NumFaces_MatchesFormula)
{
    // numFaces = (Nx+1)*Ny + Nx*(Ny+1)
    Mesh mesh;
    mesh.load(4, 4, 1.0, 1.0);
    EXPECT_EQ(mesh.numFaces(), 5 * 4 + 4 * 5);  // = 40
}

TEST(Mesh, CellIndex_RowMajorLayout)
{
    Mesh mesh;
    mesh.load(4, 8, 1.0, 1.0);
    EXPECT_EQ(mesh.cellIndex(2, 3), 2 * 8 + 3);  // = 19
}

// ── Cell centres ──────────────────────────────────────────────────────────────

TEST(Mesh, GetCellCenter_OriginCell_HalfCellOffset)
{
    // For 4×4 on [0,1]×[0,1]: dx=dy=0.25, centre(0,0)=(0.125, 0.125)
    Mesh mesh = make4x4();
    Eigen::Vector2d c = mesh.getCellCenter(0, 0);
    EXPECT_NEAR(c.x(), 0.125, 1e-12);
    EXPECT_NEAR(c.y(), 0.125, 1e-12);
}

TEST(Mesh, GetCellCenter_InternalCell_CorrectCoords)
{
    // centre(2,3) = (2.5*0.25, 3.5*0.25) = (0.625, 0.875)
    Mesh mesh = make4x4();
    Eigen::Vector2d c = mesh.getCellCenter(2, 3);
    EXPECT_NEAR(c.x(), 0.625, 1e-12);
    EXPECT_NEAR(c.y(), 0.875, 1e-12);
}

TEST(Mesh, GetCellCenter_OutOfRange_Throws)
{
    Mesh mesh = make4x4();
    EXPECT_THROW(static_cast<void>(mesh.getCellCenter(-1, 0)), std::out_of_range);
    EXPECT_THROW(static_cast<void>(mesh.getCellCenter(4,  0)), std::out_of_range);
    EXPECT_THROW(static_cast<void>(mesh.getCellCenter(0, -1)), std::out_of_range);
    EXPECT_THROW(static_cast<void>(mesh.getCellCenter(0,  4)), std::out_of_range);
}

// ── Cell volumes ──────────────────────────────────────────────────────────────

TEST(Mesh, GetCellVolume_UniformGrid_EqualsDxDy)
{
    // For 4×4 on [0,1]×[0,1]: volume = 0.25*0.25 = 0.0625
    Mesh mesh = make4x4();
    EXPECT_NEAR(mesh.getCellVolume(0),  0.0625, 1e-12);
    EXPECT_NEAR(mesh.getCellVolume(15), 0.0625, 1e-12);
}

TEST(Mesh, GetCellVolume_OutOfRange_Throws)
{
    Mesh mesh = make4x4();
    EXPECT_THROW(static_cast<void>(mesh.getCellVolume(-1)), std::out_of_range);
    EXPECT_THROW(static_cast<void>(mesh.getCellVolume(16)), std::out_of_range);
}

// ── Face areas ────────────────────────────────────────────────────────────────
// 4×4 mesh: dx=dy=0.25
// x-face area = dy = 0.25; y-face area = dx = 0.25
// xFaceCount = 5*4 = 20; first y-face at id 20.

TEST(Mesh, GetFaceArea_XFace_EqualsDy)
{
    Mesh mesh = make4x4();
    // x-face 0: i=0,j=0 → faceId=0
    EXPECT_NEAR(mesh.getFaceArea(0), 0.25, 1e-12);
    // x-face interior: i=2,j=1 → faceId=2*4+1=9
    EXPECT_NEAR(mesh.getFaceArea(9), 0.25, 1e-12);
}

TEST(Mesh, GetFaceArea_YFace_EqualsDx)
{
    Mesh mesh = make4x4();
    // First y-face: xFaceCount=20, faceId=20
    EXPECT_NEAR(mesh.getFaceArea(20), 0.25, 1e-12);
}

TEST(Mesh, GetFaceArea_OutOfRange_Throws)
{
    Mesh mesh = make4x4();
    EXPECT_THROW(static_cast<void>(mesh.getFaceArea(-1)), std::out_of_range);
    EXPECT_THROW(static_cast<void>(mesh.getFaceArea(40)), std::out_of_range);
}

// ── Face normals ──────────────────────────────────────────────────────────────

TEST(Mesh, GetFaceNormal_LeftBoundary_NegativeX)
{
    // i=0, j=0: left boundary x-face → normal=(-1,0)
    Mesh mesh = make4x4();
    Eigen::Vector2d n = mesh.getFaceNormal(0);
    EXPECT_NEAR(n.x(), -1.0, 1e-12);
    EXPECT_NEAR(n.y(),  0.0, 1e-12);
}

TEST(Mesh, GetFaceNormal_RightBoundary_PositiveX)
{
    // i=4, j=0: right boundary x-face → faceId=4*4+0=16, normal=(+1,0)
    Mesh mesh = make4x4();
    Eigen::Vector2d n = mesh.getFaceNormal(16);
    EXPECT_NEAR(n.x(), +1.0, 1e-12);
    EXPECT_NEAR(n.y(),  0.0, 1e-12);
}

TEST(Mesh, GetFaceNormal_InteriorXFace_PositiveX)
{
    // i=2, j=1 interior x-face → faceId=2*4+1=9, normal=(+1,0)
    Mesh mesh = make4x4();
    Eigen::Vector2d n = mesh.getFaceNormal(9);
    EXPECT_NEAR(n.x(), +1.0, 1e-12);
    EXPECT_NEAR(n.y(),  0.0, 1e-12);
}

TEST(Mesh, GetFaceNormal_BottomBoundary_NegativeY)
{
    // i=0, j=0 bottom y-face → faceId=xFaceCount+0=20, normal=(0,-1)
    Mesh mesh = make4x4();
    Eigen::Vector2d n = mesh.getFaceNormal(20);
    EXPECT_NEAR(n.x(),  0.0, 1e-12);
    EXPECT_NEAR(n.y(), -1.0, 1e-12);
}

TEST(Mesh, GetFaceNormal_TopBoundary_PositiveY)
{
    // i=0, j=4 top y-face → faceId=xFaceCount+0*(4+1)+4=20+4=24, normal=(0,+1)
    Mesh mesh = make4x4();
    Eigen::Vector2d n = mesh.getFaceNormal(24);
    EXPECT_NEAR(n.x(), 0.0, 1e-12);
    EXPECT_NEAR(n.y(), 1.0, 1e-12);
}

TEST(Mesh, GetFaceNormal_InteriorYFace_PositiveY)
{
    // i=1, j=2 interior y-face → faceId=20+1*(4+1)+2=27, normal=(0,+1)
    Mesh mesh = make4x4();
    Eigen::Vector2d n = mesh.getFaceNormal(27);
    EXPECT_NEAR(n.x(), 0.0, 1e-12);
    EXPECT_NEAR(n.y(), 1.0, 1e-12);
}

TEST(Mesh, GetFaceNormal_OutOfRange_Throws)
{
    Mesh mesh = make4x4();
    EXPECT_THROW(static_cast<void>(mesh.getFaceNormal(-1)), std::out_of_range);
    EXPECT_THROW(static_cast<void>(mesh.getFaceNormal(40)), std::out_of_range);
}

// ── Face neighbours ───────────────────────────────────────────────────────────

TEST(Mesh, GetNeighbors_LeftBoundary_OwnerAndNegOne)
{
    // i=0, j=0 left boundary x-face → {cell(0,0)=0, -1}
    Mesh mesh = make4x4();
    auto nb = mesh.getNeighbors(0);
    EXPECT_EQ(nb[0], 0);
    EXPECT_EQ(nb[1], -1);
}

TEST(Mesh, GetNeighbors_RightBoundary_OwnerAndNegOne)
{
    // i=4, j=2 right boundary x-face → faceId=4*4+2=18, {cell(3,2)=3*4+2=14, -1}
    Mesh mesh = make4x4();
    auto nb = mesh.getNeighbors(18);
    EXPECT_EQ(nb[0], 14);
    EXPECT_EQ(nb[1], -1);
}

TEST(Mesh, GetNeighbors_InteriorXFace_OwnerNeighbour)
{
    // i=2, j=1 interior x-face → faceId=9, owner=cell(1,1)=5, neighbour=cell(2,1)=9
    Mesh mesh = make4x4();
    auto nb = mesh.getNeighbors(9);
    EXPECT_EQ(nb[0], 5);   // cell(1,1) = 1*4+1
    EXPECT_EQ(nb[1], 9);   // cell(2,1) = 2*4+1
}

TEST(Mesh, GetNeighbors_BottomBoundary_OwnerAndNegOne)
{
    // i=0, j=0 bottom y-face → faceId=20, {cell(0,0)=0, -1}
    Mesh mesh = make4x4();
    auto nb = mesh.getNeighbors(20);
    EXPECT_EQ(nb[0], 0);
    EXPECT_EQ(nb[1], -1);
}

TEST(Mesh, GetNeighbors_TopBoundary_OwnerAndNegOne)
{
    // i=2, j=4 top y-face → faceId=20+2*(4+1)+4=34, {cell(2,3)=2*4+3=11, -1}
    Mesh mesh = make4x4();
    auto nb = mesh.getNeighbors(34);
    EXPECT_EQ(nb[0], 11);
    EXPECT_EQ(nb[1], -1);
}

TEST(Mesh, GetNeighbors_InteriorYFace_OwnerNeighbour)
{
    // i=1, j=2 interior y-face → faceId=20+1*(4+1)+2=27,
    // owner=cell(1,1)=1*4+1=5, neighbour=cell(1,2)=1*4+2=6
    Mesh mesh = make4x4();
    auto nb = mesh.getNeighbors(27);
    EXPECT_EQ(nb[0], 5);
    EXPECT_EQ(nb[1], 6);
}

TEST(Mesh, GetNeighbors_OutOfRange_Throws)
{
    Mesh mesh = make4x4();
    EXPECT_THROW(static_cast<void>(mesh.getNeighbors(-1)), std::out_of_range);
    EXPECT_THROW(static_cast<void>(mesh.getNeighbors(40)), std::out_of_range);
}
