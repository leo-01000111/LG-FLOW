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

/**
 * @brief Applies boundary conditions to a velocity or pressure field.
 *
 * Boundary conditions are applied after each field update step,
 * consistent with the SIMPLE algorithm's ghost-cell or face-value approach.
 *
 * Named patches (e.g. "top", "bottom", "left", "right") map to a
 * BoundaryPatch, so the config file drives the BC setup.
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
     * @brief Applies velocity boundary conditions to the velocity field.
     *
     * Iterates over all registered patches and enforces the prescribed
     * velocity at boundary faces using ghost-cell mirroring or direct
     * assignment depending on the boundary type.
     *
     * @param velocityField Velocity field to modify in-place.
     * @param mesh          Mesh providing face/cell topology.
     */
    void applyVelocity(Field<Eigen::Vector2d>& velocityField, const Mesh& mesh) const;

    /**
     * @brief Applies pressure boundary conditions to the pressure field.
     *
     * OUTLET faces use a fixed reference pressure (zero).
     * All other faces use zero-gradient (Neumann) conditions.
     *
     * @param pressureField Pressure field to modify in-place.
     * @param mesh          Mesh providing face/cell topology.
     */
    void applyPressure(Field<double>& pressureField, const Mesh& mesh) const;

    /**
     * @brief Returns the cell IDs of all boundary cells belonging to patches of a given type.
     *
     * Used by NavierStokesSolver to identify which cells should serve as
     * Dirichlet reference cells in the pressure-correction Poisson system
     * (e.g. all OUTLET cells pin p' = 0 so the pressure gradient can develop).
     *
     * @param type Boundary type to collect (e.g. BoundaryType::OUTLET).
     * @param mesh Mesh providing grid dimensions.
     * @return Sorted, deduplicated vector of cell IDs. Empty if no patch of the
     *         requested type is registered.
     */
    [[nodiscard]] std::vector<int> collectCellsOfType(BoundaryType type,
                                                      const Mesh&  mesh) const;

private:
    std::unordered_map<std::string, BoundaryPatch> m_patches;
};
