#pragma once

#include "core/BoundaryCondition.hpp"
#include "core/Field.hpp"
#include "core/Mesh.hpp"
#include "solver/PressureSolver.hpp"
#include "utils/Config.hpp"

#include <Eigen/Dense>

/**
 * @brief Top-level orchestrator for the incompressible Navier-Stokes solver.
 *
 * Implements the SIMPLE (Semi-Implicit Method for Pressure-Linked Equations)
 * algorithm for steady-state incompressible viscous flow:
 *
 *   Loop until convergence:
 *     1. Momentum predictor  — solve for u* using explicit pressure gradient
 *     2. Pressure correction — solve Poisson equation for p'
 *     3. Velocity correction — u = u* − (dt/ρ) ∇p'
 *     4. Apply boundary conditions
 *     5. Check residuals
 *
 * Reference: Patankar (1980), Ferziger & Perić (2020) Ch. 7.
 */
class NavierStokesSolver
{
public:
    /**
     * @brief Constructs the solver, reading all parameters from config.
     *
     * Required config keys:
     *   mesh.Nx, mesh.Ny, mesh.Lx, mesh.Ly
     *   solver.dt, solver.max_iter, solver.tolerance
     *   solver.rho, solver.nu
     *   bc.top.type, bc.top.value, bc.bottom.type, bc.left.type, bc.right.type
     *
     * @param config Parsed configuration object.
     * @throws std::runtime_error if any required key is missing.
     */
    explicit NavierStokesSolver(const Config& config);

    /**
     * @brief Allocates fields and sets initial conditions.
     *
     * Must be called before step() or run().
     * Initial condition: zero velocity everywhere, zero pressure.
     */
    void initialize();

    /**
     * @brief Advances the solution by one SIMPLE iteration.
     *
     * SIMPLE loop order (must not be reordered without explanation):
     *   1. Momentum predictor
     *   2. Pressure correction (via PressureSolver)
     *   3. Velocity correction
     *   4. Apply boundary conditions
     *
     * @param dt Time step [s] (used in the momentum predictor).
     */
    void step(double dt);

    /**
     * @brief Runs the solver until convergence or maxIter is reached.
     *
     * Calls step() in a loop, logging the residual every 100 iterations.
     * Stops when residual() < tolerance (from config).
     *
     * @param maxIter Maximum number of SIMPLE iterations.
     */
    void run(int maxIter);

    /**
     * @brief Returns the current normalised velocity residual.
     *
     * Residual = ||u_new − u_old|| / ||u_inlet|| (L2 norm, dimensionless).
     *
     * @return Normalised residual (0 = converged).
     */
    [[nodiscard]] double residual() const;

    /** @brief Read-only access to the pressure field. */
    [[nodiscard]] const Field<double>& pressure() const;

    /** @brief Read-only access to the velocity field. */
    [[nodiscard]] const Field<Eigen::Vector2d>& velocity() const;

private:
    Mesh               m_mesh;
    Field<double>      m_pressure;
    Field<Eigen::Vector2d> m_velocity;
    BoundaryCondition  m_bc;
    PressureSolver     m_pressureSolver;

    double m_dt{0.0};
    double m_rho{1.0};
    double m_nu{1e-3};
    double m_tolerance{1e-6};
    double m_residual{1.0};
};
