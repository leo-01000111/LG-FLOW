#pragma once

#include "core/Field.hpp"
#include "core/Mesh.hpp"

#include <Eigen/Dense>

/**
 * @brief Static utility class providing FVM discretization operators.
 *
 * All operators follow the cell-centred FVM convention.
 * Face fluxes are reconstructed locally and summed over each cell.
 *
 * References:
 *   Ferziger, Perić & Street — "Computational Methods for Fluid Dynamics",
 *   4th ed. (2020). Operator definitions in Chapters 4–6.
 */
class Discretization
{
public:
    Discretization() = delete;  ///< Static-only class; not instantiable.

    /**
     * @brief Computes the cell-centred divergence of a vector field.
     *
     * Discretisation: Gauss theorem, ∇·u ≈ (1/V) Σ_f (u_f · n_f) A_f
     * Face velocity u_f is obtained by linear interpolation between
     * neighbouring cell centres.
     * Scheme: central differencing (2nd order).
     *
     * @param field Vector field (e.g. velocity) [m/s].
     * @param mesh  Mesh defining cells and faces.
     * @return Scalar field of divergence values [1/s].
     */
    [[nodiscard]] static Field<double>
    divergence(const Field<Eigen::Vector2d>& field, const Mesh& mesh);

    /**
     * @brief Computes the cell-centred gradient of a scalar field.
     *
     * Discretisation: Gauss theorem, ∇φ ≈ (1/V) Σ_f φ_f n_f A_f
     * Face value φ_f obtained by linear interpolation.
     * Scheme: central differencing (2nd order).
     *
     * @param field Scalar field (e.g. pressure) [Pa].
     * @param mesh  Mesh defining cells and faces.
     * @return Vector field of gradient values [Pa/m].
     */
    [[nodiscard]] static Field<Eigen::Vector2d>
    gradient(const Field<double>& field, const Mesh& mesh);

    /**
     * @brief Computes the cell-centred Laplacian of a scalar field.
     *
     * Discretisation: ∇²φ ≈ (1/V) Σ_f (∇φ_f · n_f) A_f
     * Normal gradient at face f computed by compact two-point stencil.
     * Scheme: central differencing (2nd order).
     *
     * @param field Scalar field [Pa or m/s component].
     * @param mesh  Mesh defining cells and faces.
     * @return Scalar Laplacian field [units/m²].
     */
    [[nodiscard]] static Field<double>
    laplacian(const Field<double>& field, const Mesh& mesh);
};
