#pragma once

#include "core/BoundaryCondition.hpp"
#include "core/Field.hpp"
#include "core/Mesh.hpp"

#include <Eigen/Dense>

/**
 * @brief Computes minimum wall distance for each cell via BFS from WALL boundary cells.
 *
 *  Static utility class. Uses BFS propagation: wall-adjacent cells have distance
 *  equal to half the cell size; each subsequent layer adds the Euclidean distance
 *  between cell centres. On a uniform structured mesh this equals the exact
 *  Euclidean distance to the nearest wall face.
 *
 *  Reference: Ferziger, Peric & Street (2020) §10.1.
 */
class WallDistance
{
public:
    WallDistance() = delete;

    /**
     * @brief Computes wall distance field.
     * @param mesh  Structured mesh.
     * @param bc    Boundary condition (used to identify WALL patches).
     * @return      Field<double> with wall distance [m] at each cell centre.
     */
    [[nodiscard]] static Field<double>
    compute(const Mesh& mesh, const BoundaryCondition& bc);
};
