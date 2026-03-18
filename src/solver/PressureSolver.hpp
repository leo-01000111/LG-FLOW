#pragma once

#include "core/Field.hpp"
#include "core/Mesh.hpp"

#include <Eigen/Dense>

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

private:
    const Mesh* m_mesh;
};
