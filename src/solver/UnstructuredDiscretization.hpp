#pragma once

#include "mesh/UnstructuredField.hpp"
#include "mesh/UnstructuredMesh.hpp"

#include <Eigen/Dense>

/**
 * @brief Static FVM discretization operators for unstructured 2D meshes.
 *
 * All operators apply the Gauss theorem on the unstructured face list.
 * Interior faces use linear (central) interpolation; boundary faces use
 * the owner cell value (zero-gradient closure).
 *
 * References:
 *   Ferziger, Perić & Street — "Computational Methods for Fluid Dynamics",
 *   4th ed. (2020), Chapters 4–6.
 *   Rhie & Chow (1983), AIAA J. 21(11):1525–1532.
 */
class UnstructuredDiscretization
{
public:
    UnstructuredDiscretization() = delete;  ///< Static-only class; not instantiable.

    /**
     * @brief Gauss divergence of a vector field on an unstructured mesh.
     *
     * ∇·u ≈ (1/V) Σ_f (u_f · n_f) A_f
     * Interior face: u_f = 0.5*(u_P + u_N). Boundary face: u_f = u_owner.
     * Scheme: central differencing (2nd order) on interior faces.
     *
     * @param field Vector field (e.g. velocity) [m/s].
     * @param mesh  Unstructured mesh.
     * @return Scalar divergence field [1/s].
     */
    [[nodiscard]] static UnstructuredField<double>
    divergence(const UnstructuredField<Eigen::Vector2d>& field,
               const UnstructuredMesh& mesh);

    /**
     * @brief Gauss gradient of a scalar field on an unstructured mesh.
     *
     * ∇φ ≈ (1/V) Σ_f φ_f n_f A_f
     * Interior face: φ_f = 0.5*(φ_P + φ_N). Boundary face: φ_f = φ_owner.
     * Scheme: central differencing (2nd order) on interior faces.
     *
     * @param field Scalar field (e.g. pressure) [Pa].
     * @param mesh  Unstructured mesh.
     * @return Vector gradient field [Pa/m or m/s/m].
     */
    [[nodiscard]] static UnstructuredField<Eigen::Vector2d>
    gradient(const UnstructuredField<double>& field,
             const UnstructuredMesh& mesh);

    /**
     * @brief Rhie-Chow corrected divergence for SIMPLE pressure-velocity coupling.
     *
     * For each interior face f between owner P and neighbour N:
     *   û_f·n̂ = ū_f·n̂ − (dt/ρ)[(p_N−p_P)/|x_N−x_P| − 0.5*(∇p_P+∇p_N)·n̂]
     *
     * The correction damps checkerboard pressure modes on collocated meshes.
     * Boundary faces: owner velocity (no correction).
     *
     * Reference: Rhie & Chow (1983); Ferziger, Perić & Street (2020) Section 7.5.
     *
     * @param velocity Predicted velocity field u* [m/s].
     * @param pressure Current pressure field p [Pa].
     * @param mesh     Unstructured mesh.
     * @param dt       Time step [s].
     * @param rho      Fluid density [kg/m³].
     * @return Scalar Rhie-Chow corrected divergence [1/s].
     */
    [[nodiscard]] static UnstructuredField<double>
    divergenceRhieChow(const UnstructuredField<Eigen::Vector2d>& velocity,
                       const UnstructuredField<double>&          pressure,
                       const UnstructuredMesh&                   mesh,
                       double                                    dt,
                       double                                    rho);
};
