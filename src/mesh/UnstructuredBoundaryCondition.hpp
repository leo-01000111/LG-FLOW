#pragma once

#include "core/BoundaryCondition.hpp"
#include "mesh/UnstructuredField.hpp"
#include "mesh/UnstructuredMesh.hpp"

#include <Eigen/Dense>
#include <string>
#include <unordered_map>
#include <vector>

/**
 * @brief Applies boundary conditions to unstructured mesh fields.
 *
 * Patch names match those registered in the UnstructuredMesh (loaded from Gmsh).
 * Config keys: bc.<patchname>.type, bc.<patchname>.value_x, bc.<patchname>.value_y.
 *
 * Velocity policy:
 *   INLET / WALL      → Dirichlet: set owner cell value to patch.value.
 *   OUTLET / SYMMETRY → zero-gradient: leave cell unchanged (applied implicitly).
 *
 * Pressure policy:
 *   OUTLET → fixed reference (zero) at boundary-face owner cells.
 *   Others → zero-gradient (unchanged).
 */
class UnstructuredBoundaryCondition
{
public:
    UnstructuredBoundaryCondition() = default;

    /**
     * @brief Registers a patch with its boundary condition data.
     * @param name  Patch name matching UnstructuredMesh patch names.
     * @param patch Boundary type and prescribed value.
     */
    void addPatch(const std::string& name, const BoundaryPatch& patch);

    /**
     * @brief Returns the boundary patch for a named patch.
     * @param name Patch name.
     * @throws std::out_of_range if name is not registered.
     */
    [[nodiscard]] const BoundaryPatch& getPatch(const std::string& name) const;

    /**
     * @brief Returns true if the named patch is registered.
     * @param name Patch name.
     */
    [[nodiscard]] bool hasPatch(const std::string& name) const;

    /**
     * @brief Applies velocity boundary conditions to the velocity field.
     *
     * For each INLET/WALL patch face, sets vel[owner] = patch.value.
     *
     * @param vel  Velocity field to modify in-place.
     * @param mesh Unstructured mesh with patch topology.
     */
    void applyVelocity(UnstructuredField<Eigen::Vector2d>& vel,
                       const UnstructuredMesh& mesh) const;

    /**
     * @brief Applies pressure boundary conditions to the pressure field.
     *
     * For each OUTLET patch face, sets p[owner] = 0.0.
     *
     * @param p    Pressure field to modify in-place.
     * @param mesh Unstructured mesh with patch topology.
     */
    void applyPressure(UnstructuredField<double>& p,
                       const UnstructuredMesh& mesh) const;

    /**
     * @brief Returns CellIds owning boundary faces of the given type.
     *
     * Used to identify OUTLET cells for pressure pinning in the Poisson system.
     *
     * @param type Boundary type to collect.
     * @param mesh Mesh providing patch topology.
     * @return Sorted, deduplicated vector of owner CellIds.
     */
    [[nodiscard]] std::vector<CellId>
    collectCellsOfType(BoundaryType type, const UnstructuredMesh& mesh) const;

private:
    std::unordered_map<std::string, BoundaryPatch> m_patches;
};
