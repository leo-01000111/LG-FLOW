#pragma once

#include "core/BoundaryCondition.hpp"
#include "core/Field.hpp"
#include "core/Mesh.hpp"

#include <Eigen/Dense>
#include <vector>

/**
 * @brief Solves the pressure-correction equation in the SIMPLE algorithm.
 *
 * The SIMPLE (Semi-Implicit Method for Pressure-Linked Equations) pressure
 * correction step constructs and solves a Poisson-like equation:
 *
 *   ∇²p' = (ρ/dt) ∇·u*
 *
 * where u* is the intermediate (predicted) velocity and p' is the pressure
 * correction. After solving, pressure and velocity are corrected:
 *
 *   p  ← p  + αp  · p'
 *   u  ← u* − (dt/ρ) ∇p'
 *
 * Reference: Patankar (1980) "Numerical Heat Transfer and Fluid Flow", Ch. 6.
 */
class PressureSolver
{
public:
    /**
     * @brief Constructs the pressure solver bound to a mesh.
     * @param mesh Reference to the computational mesh (must outlive this object).
     */
    explicit PressureSolver(const Mesh& mesh);

    /**
     * @brief Solves for the pressure correction and updates pressure and velocity.
     *
     * Steps performed:
     *   1. Assemble the pressure-correction matrix from ∇·u*.
     *   2. Solve the linear system (Gauss-Seidel or conjugate gradient).
     *   3. Update pressure:  p  += αp * p'
     *   4. Update velocity:  u  -= (dt/ρ) * ∇p'
     *
     * @param velocityField  Predicted velocity u* (modified in-place).
     * @param pressureField  Current pressure (modified in-place).
     * @param dt             Time step [s].
     * @param rho            Fluid density [kg/m³].
     * @param alphaP         Pressure under-relaxation factor (0 < αp ≤ 1).
     * @return               L2 norm of the pressure-correction residual [Pa].
     */
    [[nodiscard]] double solve(Field<Eigen::Vector2d>& velocityField,
                               Field<double>&          pressureField,
                               double                  dt,
                               double                  rho,
                               double                  alphaP);

    /**
     * @brief Supplies the boundary closure used by every subsequent solve().
     *
     * Without it, every boundary face is zero-gradient and cell 0 is pinned
     * (legacy closure, kept for standalone use and existing tests). When set:
     *   - The source uses the boundary-aware Rhie-Chow divergence (exact
     *     wall/inlet fluxes, Rhie-Chow outlet faces).
     *   - FIXED pressure faces (outlets) impose p' = 0 on the face, adding
     *     A_f / d_b to the owner's diagonal (d_b = half-cell distance).
     *   - Without any FIXED face (all-wall cavity) cell 0 is pinned to remove
     *     the pure-Neumann singularity.
     *   - The velocity correction uses a Gauss gradient of p' with p' = 0 on
     *     FIXED faces and zero-gradient elsewhere.
     * Reference: Ferziger, Perić & Street (2020) §7.5 and Ch. 8.
     *
     * @param pressureFaces Pressure closure for every boundary face.
     * @param velocityFaces Velocity closure, same face order as pressureFaces.
     * @throws std::invalid_argument if the lists differ in size.
     */
    void setBoundaryFaces(std::vector<PressureBoundaryFace> pressureFaces,
                          std::vector<VelocityBoundaryFace> velocityFaces);

private:
    const Mesh* m_mesh;

    std::vector<PressureBoundaryFace> m_pressureFaces;    ///< Empty → legacy closure
    std::vector<VelocityBoundaryFace> m_velocityFaces;    ///< Same order as m_pressureFaces
    std::vector<PressureBoundaryFace> m_correctionFaces;  ///< m_pressureFaces with p'_b = 0
    bool m_hasFixedFace{false};  ///< At least one FIXED face → no reference cell needed
};
