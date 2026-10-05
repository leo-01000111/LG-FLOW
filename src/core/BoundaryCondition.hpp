#pragma once

#include "core/Field.hpp"
#include "core/Mesh.hpp"

#include <Eigen/Dense>
#include <string>
#include <unordered_map>
#include <vector>

/**
 * @brief Identifies which physical boundary a face group belongs to.
 */
enum class BoundaryType
{
    INLET,           ///< Prescribed uniform velocity (Dirichlet)
    PARABOLIC_INLET, ///< Poiseuille profile: u(y)=6·U_avg·y·(Ly-y)/Ly²; value.x() = U_avg
    OUTLET,          ///< Zero-gradient velocity, fixed (zero) pressure
    WALL,            ///< No-slip: zero velocity at the face
    SYMMETRY         ///< Zero normal gradient for all quantities
};

/**
 * @brief Stores the type and prescribed value for one boundary patch.
 *
 * For INLET: m_value is the inlet velocity [m/s].
 * For WALL:  m_value is the wall velocity (zero for stationary, non-zero for lid).
 * For OUTLET/SYMMETRY: m_value is unused.
 */
struct BoundaryPatch
{
    BoundaryType     type;
    Eigen::Vector2d  value{Eigen::Vector2d::Zero()};  ///< Prescribed vector value
};

/** @brief How the velocity is closed at a boundary face. */
enum class VelocityBCKind
{
    DIRICHLET,     ///< Face velocity prescribed (WALL, INLET, PARABOLIC_INLET)
    ZERO_GRADIENT  ///< Face velocity = owner-cell velocity (OUTLET, SYMMETRY, unpatched)
};

/** @brief How the pressure is closed at a boundary face. */
enum class PressureBCKind
{
    FIXED,         ///< Face pressure prescribed (OUTLET, p = 0)
    ZERO_GRADIENT  ///< Face pressure = owner-cell pressure (all other types)
};

/**
 * @brief Precomputed velocity boundary data for one boundary face.
 *
 * Built once by BoundaryCondition::buildVelocityFaces() and handed to the
 * consumers that need it (MomentumSolver, Rhie-Chow divergence), so the hot
 * loops never look up patches by name.
 */
struct VelocityBoundaryFace
{
    int             face{-1};                          ///< Global face ID
    int             owner{-1};                         ///< Interior cell owning the face
    Eigen::Vector2d normal{Eigen::Vector2d::Zero()};   ///< Outward unit normal
    double          area{0.0};                         ///< Face area (length in 2D) [m]
    double          distance{0.0};                     ///< Owner centre → face centre [m]
    VelocityBCKind  kind{VelocityBCKind::ZERO_GRADIENT};
    Eigen::Vector2d value{Eigen::Vector2d::Zero()};    ///< Face velocity if DIRICHLET [m/s]
};

/**
 * @brief Precomputed pressure boundary data for one boundary face.
 *
 * Built once by BoundaryCondition::buildPressureFaces() and handed to the
 * consumers that need it (PressureSolver, pressure gradient, Rhie-Chow).
 */
struct PressureBoundaryFace
{
    int             face{-1};                          ///< Global face ID
    int             owner{-1};                         ///< Interior cell owning the face
    Eigen::Vector2d normal{Eigen::Vector2d::Zero()};   ///< Outward unit normal
    double          area{0.0};                         ///< Face area (length in 2D) [m]
    double          distance{0.0};                     ///< Owner centre → face centre [m]
    PressureBCKind  kind{PressureBCKind::ZERO_GRADIENT};
    double          value{0.0};                        ///< Face pressure if FIXED [Pa]
};

/**
 * @brief Holds the boundary patches and turns them into boundary-face closures.
 *
 * Named patches (e.g. "top", "bottom", "left", "right") map to a
 * BoundaryPatch, so the config file drives the BC setup.
 *
 * The solver applies boundary conditions as boundary-face fluxes:
 * buildVelocityFaces() / buildPressureFaces() produce precomputed per-face
 * data that NavierStokesSolver hands to the momentum and pressure solvers once
 * in initialize(). The cell-overwrite methods applyVelocity()/applyPressure()
 * are legacy and no longer used by NavierStokesSolver.
 */
class BoundaryCondition
{
public:
    BoundaryCondition() = default;

    /**
     * @brief Registers a named boundary patch.
     * @param name  Patch name matching the config key (e.g. "top").
     * @param patch The boundary type and prescribed value.
     */
    void addPatch(const std::string& name, const BoundaryPatch& patch);

    /**
     * @brief Returns the patch data for a named boundary.
     * @param name Patch name.
     * @return Const reference to the BoundaryPatch.
     * @throws std::out_of_range if name is not registered.
     */
    [[nodiscard]] const BoundaryPatch& getPatch(const std::string& name) const;

    /**
     * @brief Legacy: overwrites boundary cells with the prescribed velocity.
     *
     * Sets the boundary *cells* (not faces) of WALL / INLET / PARABOLIC_INLET
     * patches. This places walls at the first cell centres instead of the
     * domain edge, so NavierStokesSolver no longer calls it; boundary
     * conditions enter through buildVelocityFaces() instead.
     *
     * @param velocityField Velocity field to modify in-place.
     * @param mesh          Mesh providing face/cell topology.
     */
    void applyVelocity(Field<Eigen::Vector2d>& velocityField, const Mesh& mesh) const;

    /**
     * @brief Legacy: sets OUTLET boundary cells to zero pressure.
     *
     * Superseded by buildPressureFaces() (p = 0 on the outlet face); no longer
     * called by NavierStokesSolver.
     *
     * @param pressureField Pressure field to modify in-place.
     * @param mesh          Mesh providing face/cell topology.
     */
    void applyPressure(Field<double>& pressureField, const Mesh& mesh) const;

    /**
     * @brief Returns the cell IDs of all boundary cells belonging to patches of a given type.
     *
     * COMPATIBILITY PATH (Spalart-Allmaras, WallDistance): "wall cells" are
     * the cells that own a WALL face. WallDistance measures distances to their
     * centres (offset by half a cell from the true wall at the face) and the
     * SA model imposes ν̃ = 0 on them rather than on the wall faces. Moving
     * these to the faces is deferred work.
     *
     * @param type Boundary type to collect (e.g. BoundaryType::OUTLET).
     * @param mesh Mesh providing grid dimensions.
     * @return Sorted, deduplicated vector of cell IDs. Empty if no patch of the
     *         requested type is registered.
     */
    [[nodiscard]] std::vector<int> collectCellsOfType(BoundaryType type,
                                                      const Mesh&  mesh) const;

    /**
     * @brief Builds the velocity closure for every boundary face of the mesh.
     *
     * Mapping: WALL / INLET → DIRICHLET with the patch value; PARABOLIC_INLET →
     * DIRICHLET with u(y) = 6·U_avg·y·(Ly−y)/Ly² evaluated at the face centre;
     * OUTLET, SYMMETRY and sides without a patch → ZERO_GRADIENT.
     *
     * Entries are in ascending global face order, identical to
     * buildPressureFaces(), so index k refers to the same face in both.
     *
     * @param mesh Mesh providing face topology and geometry.
     * @return One entry per boundary face.
     * @throws std::invalid_argument if PARABOLIC_INLET is used on "top" or "bottom".
     */
    [[nodiscard]] std::vector<VelocityBoundaryFace> buildVelocityFaces(const Mesh& mesh) const;

    /**
     * @brief Builds the pressure closure for every boundary face of the mesh.
     *
     * Mapping: OUTLET → FIXED with p = 0; all other types and unpatched sides →
     * ZERO_GRADIENT. Same face order as buildVelocityFaces().
     *
     * @param mesh Mesh providing face topology and geometry.
     * @return One entry per boundary face.
     */
    [[nodiscard]] std::vector<PressureBoundaryFace> buildPressureFaces(const Mesh& mesh) const;

private:
    std::unordered_map<std::string, BoundaryPatch> m_patches;
};
