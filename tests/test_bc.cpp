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

// ── collectCellsOfType ────────────────────────────────────────────────────────

TEST(BoundaryCondition, CollectCellsOfType_NoPatchOfType_ReturnsEmpty)
{
    // No OUTLET patch registered → must return an empty vector.
    Mesh mesh = makeSquareMesh();  // 8x8

    BoundaryCondition bc;
    bc.addPatch("top",    {BoundaryType::WALL, Eigen::Vector2d::Zero()});
    bc.addPatch("bottom", {BoundaryType::WALL, Eigen::Vector2d::Zero()});
    bc.addPatch("left",   {BoundaryType::INLET, Eigen::Vector2d(1.0, 0.0)});

    const std::vector<int> cells = bc.collectCellsOfType(BoundaryType::OUTLET, mesh);
    EXPECT_TRUE(cells.empty());
}

TEST(BoundaryCondition, CollectCellsOfType_RightOutlet_CorrectCellCount)
{
    // An 8x8 mesh has Ny=8 cells in the right column → 8 outlet cells.
    Mesh mesh = makeSquareMesh();  // 8x8, Ny=8

    BoundaryCondition bc;
    bc.addPatch("right", {BoundaryType::OUTLET, Eigen::Vector2d::Zero()});

    const std::vector<int> cells = bc.collectCellsOfType(BoundaryType::OUTLET, mesh);

    ASSERT_EQ(static_cast<int>(cells.size()), mesh.Ny())
        << "Expected Ny outlet cells for right OUTLET patch";

    // Cell IDs for right column: (Nx-1)*Ny + j, j in [0, Ny)
    const int Nx = mesh.Nx();
    const int Ny = mesh.Ny();
    for (int j = 0; j < Ny; ++j)
    {
        const int expected = (Nx - 1) * Ny + j;
        EXPECT_NE(std::find(cells.begin(), cells.end(), expected), cells.end())
            << "Cell (Nx-1, " << j << ") id=" << expected << " not found";
    }
}

TEST(BoundaryCondition, CollectCellsOfType_ResultIsSortedAndUnique)
{
    // Collecting OUTLET for a single patch must return a sorted, deduplicated list.
    Mesh mesh = makeSquareMesh();  // 8x8

    BoundaryCondition bc;
    bc.addPatch("right", {BoundaryType::OUTLET, Eigen::Vector2d::Zero()});

    const std::vector<int> cells = bc.collectCellsOfType(BoundaryType::OUTLET, mesh);

    for (std::size_t k = 1; k < cells.size(); ++k)
        EXPECT_LT(cells[k - 1], cells[k])
            << "Result is not strictly sorted at index " << k;
}

TEST(BoundaryCondition, CollectCellsOfType_LeftInlet_CorrectCellIds)
{
    // Left column cells: i=0, j in [0, Ny) → cell ID = 0 * Ny + j = j
    Mesh mesh = makeSquareMesh();  // 8x8, Ny=8

    BoundaryCondition bc;
    bc.addPatch("left", {BoundaryType::INLET, Eigen::Vector2d(1.0, 0.0)});

    const std::vector<int> cells = bc.collectCellsOfType(BoundaryType::INLET, mesh);

    ASSERT_EQ(static_cast<int>(cells.size()), mesh.Ny());
    for (int j = 0; j < mesh.Ny(); ++j)
        EXPECT_EQ(cells[static_cast<std::size_t>(j)], j)
            << "Left inlet cell j=" << j << " should have id=" << j;
}

// ── buildVelocityFaces / buildPressureFaces ───────────────────────────────────

namespace
{
/// Channel mesh: 8×4 cells on [0, 2] × [0, 1] → dx = 0.25, dy = 0.25.
Mesh makeChannelMesh()
{
    Mesh m;
    m.load(8, 4, 2.0, 1.0);
    return m;
}

BoundaryCondition makeChannelBC(BoundaryType inletType)
{
    BoundaryCondition bc;
    bc.addPatch("left",   {inletType,            Eigen::Vector2d(1.0, 0.0)});
    bc.addPatch("right",  {BoundaryType::OUTLET, Eigen::Vector2d::Zero()});
    bc.addPatch("bottom", {BoundaryType::WALL,   Eigen::Vector2d::Zero()});
    bc.addPatch("top",    {BoundaryType::WALL,   Eigen::Vector2d::Zero()});
    return bc;
}
}  // namespace

TEST(BoundaryCondition, BuildVelocityFaces_Channel_OneEntryPerBoundaryFace)
{
    const Mesh mesh = makeChannelMesh();
    const BoundaryCondition bc = makeChannelBC(BoundaryType::INLET);

    const std::vector<VelocityBoundaryFace> faces = bc.buildVelocityFaces(mesh);

    // 2·Ny x-faces + 2·Nx y-faces = 8 + 16
    EXPECT_EQ(static_cast<int>(faces.size()), 2 * mesh.Ny() + 2 * mesh.Nx());
    for (std::size_t k = 1; k < faces.size(); ++k)
        EXPECT_LT(faces[k - 1].face, faces[k].face) << "not in ascending face order";
}

TEST(BoundaryCondition, BuildVelocityFaces_Channel_KindsMatchPatchTypes)
{
    const Mesh mesh = makeChannelMesh();
    const BoundaryCondition bc = makeChannelBC(BoundaryType::INLET);

    for (const VelocityBoundaryFace& f : bc.buildVelocityFaces(mesh))
    {
        if (f.normal.x() > 0.5)  // right: outlet
        {
            EXPECT_EQ(f.kind, VelocityBCKind::ZERO_GRADIENT);
        }
        else if (f.normal.x() < -0.5)  // left: uniform inlet
        {
            EXPECT_EQ(f.kind, VelocityBCKind::DIRICHLET);
            EXPECT_NEAR(f.value.x(), 1.0, 1e-12);
            EXPECT_NEAR(f.value.y(), 0.0, 1e-12);
        }
        else  // top / bottom: stationary walls
        {
            EXPECT_EQ(f.kind, VelocityBCKind::DIRICHLET);
            EXPECT_NEAR(f.value.norm(), 0.0, 1e-12);
        }
    }
}

TEST(BoundaryCondition, BuildVelocityFaces_Channel_DistanceIsHalfCell)
{
    const Mesh mesh = makeChannelMesh();  // dx = dy = 0.25
    const BoundaryCondition bc = makeChannelBC(BoundaryType::INLET);

    for (const VelocityBoundaryFace& f : bc.buildVelocityFaces(mesh))
    {
        EXPECT_NEAR(f.distance, 0.125, 1e-12) << "face " << f.face;
        EXPECT_NEAR(f.area,     0.25,  1e-12) << "face " << f.face;
    }
}

TEST(BoundaryCondition, BuildVelocityFaces_ParabolicInlet_FluxEqualsUavgTimesH)
{
    const Mesh mesh = makeChannelMesh();  // H = 1, U_avg = 1
    const BoundaryCondition bc = makeChannelBC(BoundaryType::PARABOLIC_INLET);

    double inflow = 0.0;
    for (const VelocityBoundaryFace& f : bc.buildVelocityFaces(mesh))
        if (f.normal.x() < -0.5)
            inflow += f.value.x() * f.area;

    // Face values are face averages of the profile, so the discrete inflow is
    // exact (a point-sampled profile would give U·H·(1 + dy²/(2H²)) = 1.03125).
    EXPECT_NEAR(inflow, 1.0, 1e-12);
}

TEST(BoundaryCondition, BuildVelocityFaces_ParabolicOnTop_Throws)
{
    const Mesh mesh = makeChannelMesh();
    BoundaryCondition bc;
    bc.addPatch("top", {BoundaryType::PARABOLIC_INLET, Eigen::Vector2d(1.0, 0.0)});

    EXPECT_THROW(static_cast<void>(bc.buildVelocityFaces(mesh)), std::invalid_argument);
}

TEST(BoundaryCondition, BuildPressureFaces_Channel_OnlyOutletIsFixed)
{
    const Mesh mesh = makeChannelMesh();
    const BoundaryCondition bc = makeChannelBC(BoundaryType::INLET);

    for (const PressureBoundaryFace& f : bc.buildPressureFaces(mesh))
    {
        if (f.normal.x() > 0.5)
        {
            EXPECT_EQ(f.kind, PressureBCKind::FIXED);
            EXPECT_NEAR(f.value, 0.0, 1e-12);
        }
        else
        {
            EXPECT_EQ(f.kind, PressureBCKind::ZERO_GRADIENT);
        }
    }
}

TEST(BoundaryCondition, BuildFaces_SameMesh_VelocityAndPressureOrderMatch)
{
    const Mesh mesh = makeChannelMesh();
    const BoundaryCondition bc = makeChannelBC(BoundaryType::INLET);

    const std::vector<VelocityBoundaryFace> vf = bc.buildVelocityFaces(mesh);
    const std::vector<PressureBoundaryFace> pf = bc.buildPressureFaces(mesh);

    ASSERT_EQ(vf.size(), pf.size());
    for (std::size_t k = 0; k < vf.size(); ++k)
    {
        EXPECT_EQ(vf[k].face,  pf[k].face);
        EXPECT_EQ(vf[k].owner, pf[k].owner);
    }
}

TEST(BoundaryCondition, BuildVelocityFaces_NoPatches_AllZeroGradient)
{
    const Mesh mesh = makeChannelMesh();
    const BoundaryCondition bc;

    for (const VelocityBoundaryFace& f : bc.buildVelocityFaces(mesh))
        EXPECT_EQ(f.kind, VelocityBCKind::ZERO_GRADIENT);
}
