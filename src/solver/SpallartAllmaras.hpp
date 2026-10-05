#pragma once

#include "core/Field.hpp"
#include "core/Mesh.hpp"

#include <Eigen/Dense>

/**
 * @brief Spalart-Allmaras one-equation turbulence model.
 *
 *  Solves the transport equation for ν̃ (modified eddy viscosity):
 *
 *    ∂ν̃/∂t + u·∇ν̃ = cb1·S̃·ν̃
 *                   - cw1·fw·(ν̃/d)²
 *                   + (1/σ)·∇·((ν+ν̃)∇ν̃)
 *                   + (cb2/σ)·(∇ν̃)²
 *
 *  Reference: Spalart, P.R. & Allmaras, S.R. (1994).
 *    "A one-equation turbulence model for aerodynamic flows."
 *    La Recherche Aérospatiale, 1, 5–21.
 *    Allmaras, S.R., Johnson, F.T. & Spalart, P.R. (2012).
 *    Modifications and clarifications for implementation of the SA model.
 *    ICCFD7-1902.
 */
class SpallartAllmaras
{
public:
    /** @brief Constructs the SA model bound to a mesh.
     *  @param mesh Mesh (must outlive this object).
     */
    explicit SpallartAllmaras(const Mesh& mesh);

    // ── Closure functions (all static, all nodiscard) ─────────────────────

    /** @brief χ = ν̃/ν */
    [[nodiscard]] static double chi(double nuTilde, double nu);

    /** @brief fv1 = χ³/(χ³ + cv1³),  cv1=7.1  (wall damping function) */
    [[nodiscard]] static double fv1(double chiVal);

    /** @brief fv2 = 1 - χ/(1 + χ·fv1)  (SA-2003 Allmaras correction) */
    [[nodiscard]] static double fv2(double chiVal, double fv1Val);

    /** @brief S̃ = Ω + (ν̃/(κ²d²))·fv2  (modified strain rate)
     *  @param omega    Vorticity magnitude [1/s]
     *  @param d        Wall distance [m]
     *  @param nu       Molecular viscosity [m²/s]
     *  @param nuTilde  SA variable [m²/s]
     */
    [[nodiscard]] static double sTilde(double omega, double d, double nu, double nuTilde);

    /** @brief r = ν̃/(S̃·κ²·d²), clamped to r_lim = 10 */
    [[nodiscard]] static double rFunc(double nuTilde, double sTildeVal, double d);

    /** @brief g = r + cw2·(r⁶ - r) */
    [[nodiscard]] static double gFunc(double rVal);

    /** @brief fw = g·[(1+cw3⁶)/(g⁶+cw3⁶)]^(1/6)  (log-layer shaping) */
    [[nodiscard]] static double fw(double gVal);

    // ── Transport equation ────────────────────────────────────────────────

    /** @brief Advances the SA transport equation by one step.
     *
     *  Builds and solves a sparse linear system (same SparseLU pattern as
     *  MomentumSolver). Discretisation:
     *    - Transient: Euler implicit
     *    - Diffusion: implicit central, coefficient (ν+ν̃)/σ per face
     *    - Convection: first-order upwind (explicit, added to RHS)
     *    - Production: explicit source  cb1·S̃·ν̃
     *    - Destruction: linearized implicit  cw1·fw·ν̃^{k}/d²  (multiplied into diagonal)
     *    - Cross-diffusion: explicit  (cb2/σ)·(∇ν̃)²
     *
     *  After solve, negative ν̃ values are clamped to 0 (Allmaras 2012 §2.4).
     *  Under-relaxation is applied: ν̃_new = α·ν̃_solved + (1-α)·ν̃_old.
     *
     *  @param nuTilde   SA variable field [m²/s] — updated in place.
     *  @param velocity  Current velocity field u^k.
     *  @param wallDist  Wall distance field d [m].
     *  @param nu        Molecular kinematic viscosity [m²/s].
     *  @param dt        Time step [s].
     *  @param alphaNu   Under-relaxation factor in (0,1].
     *  @throws std::runtime_error on non-finite ν̃ after solve.
     */
    void solve(Field<double>&                nuTilde,
               const Field<Eigen::Vector2d>& velocity,
               const Field<double>&          wallDist,
               double                        nu,
               double                        dt,
               double                        alphaNu = 0.7);

    /** @brief Computes turbulent viscosity ν_t = ν̃·fv1(χ).
     *  Negative ν̃ → ν_t = 0 (Allmaras 2012 fix, no negative branch).
     *  @param nuTilde  SA variable field.
     *  @param nu       Molecular viscosity [m²/s].
     *  @param mesh     Mesh the field lives on.
     *  @return         ν_t field [m²/s].
     */
    [[nodiscard]] static Field<double>
    computeNuT(const Field<double>& nuTilde, double nu, const Mesh& mesh);

private:
    const Mesh* m_mesh;
};
