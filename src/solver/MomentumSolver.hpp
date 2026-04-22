#pragma once

#include "core/Field.hpp"
#include "core/Mesh.hpp"

#include <Eigen/Dense>

/**
 * @brief Convection discretization options.
 *
 * Defined here so both MomentumSolver and NavierStokesSolver share the same enum.
 * UPWIND:  First-order donor-cell (upwind) — unconditionally stable.
 * CENTRAL: Second-order Gauss face-average — more accurate, less dissipative.
 */
enum class ConvectionScheme { CENTRAL, UPWIND };

/**
 * @brief Solves the semi-implicit momentum predictor step of SIMPLE.
 *
 * For each velocity component φ ∈ {ux, uy}, assembles and solves:
 *
 *   (V_P/dt + ν Σ_f a_f) φ*_P − ν Σ_f a_f φ*_nb = b_P
 *
 * where:
 *   a_f   = A_f / |x_N − x_P|  (face diffusion coefficient)
 *   b_P   = (V_P/dt) φ^k_P  −  V_P (u·∇)φ|^k_P  −  V_P (1/ρ) ∂p^k/∂φ_P
 *
 * - Viscous diffusion: implicit (central, 2nd order) — removes explicit CFL constraint.
 * - Convection: explicit (upwind or central, selected per call).
 * - Boundary faces: zero-normal-gradient closure (no off-diagonal coupling to ghost cells).
 *   Dirichlet BCs are enforced by NavierStokesSolver after the solve.
 * - Linear solver: SparseLU (direct, exact for the 5-point structured stencil).
 *
 * Reference: Ferziger, Perić & Street, 4th ed. (2020), Section 7.4 (eq. 7.17).
 */
class MomentumSolver
{
public:
    /**
     * @brief Constructs bound to a mesh.
     * @param mesh Mesh (must outlive this object).
     */
    explicit MomentumSolver(const Mesh& mesh);

    /**
     * @brief Solves the momentum predictor; writes result to uStar.
     *
     * Under-relaxation is embedded into the matrix (Patankar 1980, eq. 6.36):
     *   a_P_mod = a_P / α_u
     *   b_mod   = b + (1 − α_u)/α_u · a_P · u^k_P
     * This preserves mass conservation in u* so the pressure-corrected velocity
     * can be committed directly without a post-solve mixing step.
     *
     * @param uStar     Output: intermediate velocity u* (written in place).
     * @param uOld      Velocity field u^k from the previous SIMPLE iteration.
     * @param pressure  Pressure field p^k from the previous SIMPLE iteration.
     * @param dt        Effective time step [s]. Must be > 0.
     * @param rho       Fluid density [kg/m³]. Must be > 0.
     * @param nu        Kinematic viscosity [m²/s]. Must be >= 0.
     * @param scheme    Convection scheme (UPWIND or CENTRAL).
     * @param alphaU    Velocity under-relaxation factor in (0, 1]. Default 1.0 (no relaxation).
     * @throws std::invalid_argument on bad parameter values.
     * @throws std::runtime_error if the linear solve fails or produces non-finite values.
     */
    void solve(Field<Eigen::Vector2d>&       uStar,
               const Field<Eigen::Vector2d>& uOld,
               const Field<double>&          pressure,
               double                        dt,
               double                        rho,
               double                        nu,
               ConvectionScheme              scheme,
               double                        alphaU = 1.0);

private:
    const Mesh* m_mesh;
};
