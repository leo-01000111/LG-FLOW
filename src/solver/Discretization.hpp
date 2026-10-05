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

    /**
     * @brief Computes the cell-centred divergence with Rhie-Chow face interpolation.
     *
     * Replaces the face velocity used in the continuity equation with a
     * Rhie-Chow-corrected value to prevent pressure-velocity decoupling
     * (checkerboard oscillations) on collocated meshes.
     *
     * For each interior face f between owner P and neighbour N, the normal
     * face velocity is:
     *
     *   û_f · n̂_f = ū_f · n̂_f
     *             − (dt/ρ) [(p_N − p_P)/|x_N − x_P|
     *                       − 0.5·((∇p)_P + (∇p)_N) · n̂_f]
     *
     * where ū_f = 0.5(u_P + u_N) is the linearly interpolated face velocity.
     * The correction term damps checkerboard pressure modes by restoring
     * sensitivity to the compact (two-point) pressure gradient.
     *
     * Boundary faces: same as divergence() — owner value, no correction.
     *
     * Reference: Rhie & Chow (1983), AIAA J. 21(11):1525–1532.
     *            Ferziger, Perić & Street (2020) Section 7.5.
     *
     * @param velocity  Predicted velocity field u* [m/s].
     * @param pressure  Current pressure field p^k [Pa].
     * @param mesh      Mesh defining cells and faces.
     * @param dt        Effective time step [s].
     * @param rho       Fluid density [kg/m³].
     * @return Scalar field of Rhie-Chow-corrected divergence values [1/s].
     */
    [[nodiscard]] static Field<double>
    divergenceRhieChow(const Field<Eigen::Vector2d>& velocity,
                       const Field<double>&          pressure,
                       const Mesh&                   mesh,
                       double                        dt,
                       double                        rho);
};
