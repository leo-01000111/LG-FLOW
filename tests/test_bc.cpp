#include "core/BoundaryCondition.hpp"
#include "core/Field.hpp"
#include "core/Mesh.hpp"

#include <Eigen/Dense>
#include <gtest/gtest.h>
#include <stdexcept>

// ── Milestone 3 will replace these stubs with real assertions ────────────────

namespace
{
Mesh makeSquareMesh()
{
    Mesh m;
    m.load(8, 8, 1.0, 1.0);
    return m;
}
}  // namespace

/**
 * When implemented: addPatch + getPatch should retrieve the same BoundaryPatch.
 */
TEST(BoundaryCondition, AddAndGetPatch_RoundTrip)
{
    BoundaryCondition bc;
    BoundaryPatch patch{BoundaryType::WALL, Eigen::Vector2d(1.0, 0.0)};
    bc.addPatch("top", patch);

    const BoundaryPatch& retrieved = bc.getPatch("top");
    EXPECT_EQ(retrieved.type, BoundaryType::WALL);
    // TODO(Milestone 3): EXPECT_NEAR(retrieved.value.x(), 1.0, 1e-12)
    EXPECT_TRUE(true);
}

/**
 * When implemented: getPatch on an unknown name should throw std::out_of_range.
 */
TEST(BoundaryCondition, GetPatch_UnknownName_Throws)
{
    BoundaryCondition bc;
    EXPECT_THROW(static_cast<void>(bc.getPatch("nonexistent")), std::out_of_range);
}

/**
 * When implemented: applyVelocity with a WALL BC should set the top-boundary
 * cells to the lid velocity (1, 0) for a lid-driven cavity.
 */
TEST(BoundaryCondition, ApplyVelocity_LidWall_SetsLidVelocity)
{
    Mesh mesh = makeSquareMesh();
    Field<Eigen::Vector2d> vel(mesh, Eigen::Vector2d::Zero());

    BoundaryCondition bc;
    bc.addPatch("top",    {BoundaryType::WALL,   Eigen::Vector2d(1.0, 0.0)});
    bc.addPatch("bottom", {BoundaryType::WALL,   Eigen::Vector2d::Zero()});
    bc.addPatch("left",   {BoundaryType::WALL,   Eigen::Vector2d::Zero()});
    bc.addPatch("right",  {BoundaryType::WALL,   Eigen::Vector2d::Zero()});

    EXPECT_NO_THROW(bc.applyVelocity(vel, mesh));

    // TODO(Milestone 3): verify that top boundary cells have vel = (1,0)
    EXPECT_TRUE(true);
}

/**
 * When implemented: applyPressure with an OUTLET patch should set the
 * corresponding boundary cells to zero pressure.
 */
TEST(BoundaryCondition, ApplyPressure_Outlet_ZeroPressureBoundary)
{
    Mesh mesh = makeSquareMesh();
    Field<double> p(mesh, 1.0);  // start with non-zero pressure

    BoundaryCondition bc;
    bc.addPatch("right", {BoundaryType::OUTLET, Eigen::Vector2d::Zero()});

    EXPECT_NO_THROW(bc.applyPressure(p, mesh));

    // TODO(Milestone 3): verify that right boundary cells have p ≈ 0
    EXPECT_TRUE(true);
}
