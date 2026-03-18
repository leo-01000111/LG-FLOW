#pragma once

#include "core/BoundaryCondition.hpp"
#include "core/Field.hpp"
#include "core/Mesh.hpp"
#include "io/VTKWriter.hpp"
#include "solver/PressureSolver.hpp"
#include "utils/Config.hpp"

#include <Eigen/Dense>
#include <optional>
#include <string>

/**
 * @brief Top-level orchestrator for the incompressible Navier-Stokes solver.
 *
 * Implements the SIMPLE (Semi-Implicit Method for Pressure-Linked Equations)
 * algorithm for steady-state incompressible viscous flow:
 *
 *   Loop until convergence:
 *     1. Momentum predictor  -- solve for u* using explicit pressure gradient
 *     2. Pressure correction -- solve Poisson equation for p'
 *     3. Velocity correction -- u = u* - (dt/rho) * grad(p')
 *     4. Apply boundary conditions
 *     5. Check residuals
 *
 * Reference: Patankar (1980), Ferziger & Peric (2020) Ch. 7.
 *
 * ### Initialization contract
 *
 * The object has two distinct lifecycle phases:
 *
 *   1. **Constructed** (after constructor): config parameters are loaded;
 *      no Mesh, Field, or PressureSolver storage is allocated.
 *
 *   2. **Initialized** (after initialize()): Mesh is built, Fields and
 *      PressureSolver are allocated and zeroed. All methods are safe to call.
 *
 * Calling pressure(), velocity(), step(), or run() before initialize() throws
 * std::logic_error with an actionable message.
 */
class NavierStokesSolver
{
public:
    /**
     * @brief Constructs the solver, reading all parameters from config.
     *
     * Only config parameters are read here; no mesh or field storage is
     * allocated. Call initialize() before using the solver.
     *
     * Required config keys (all have safe defaults):
     *   mesh.Nx (default 16), mesh.Ny (default 16)
     *   mesh.Lx (default 1.0), mesh.Ly (default 1.0)
     *   solver.dt (default 0.01), solver.rho (default 1.0)
     *   solver.nu (default 0.01), solver.tolerance (default 1e-6)
     *
     * @param config Parsed configuration object.
     */
    explicit NavierStokesSolver(const Config& config);

    /**
     * @brief Allocates the mesh, fields, and pressure solver; sets initial
     *        conditions to zero.
     *
     * Must be called before step(), run(), pressure(), or velocity().
     * Calling it more than once re-initializes everything (mesh rebuild,
     * fields zeroed).
     *
     * @throws std::invalid_argument (propagated from Mesh::load()) if config
     *         mesh parameters are invalid.
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
     * @param dt Time step [s].
     * @throws std::logic_error if called before initialize().
     */
    void step(double dt);

    /**
     * @brief Runs the solver until convergence or maxIter is reached.
     *
     * Calls step() in a loop, logging progress every 50 iterations.
     * Stops when velocityResidual() < tolerance AND
     * continuityResidual() < tolerance (from config).
     * Writes history.csv in output.dir with one row per iteration.
     * Writes VTK snapshots every output.vtk_interval iterations.
     *
     * @param maxIter Maximum number of SIMPLE iterations.
     * @throws std::logic_error   if called before initialize().
     * @throws std::runtime_error if history.csv cannot be created.
     */
    void run(int maxIter);

    /**
     * @brief Returns the current normalised velocity residual.
     *
     * Equivalent to velocityResidual(). Retained for API compatibility.
     *
     * @return ||u_new - u_old|| / max(||u_old||, 1e-12)
     */
    [[nodiscard]] double residual() const;

    /**
     * @brief Returns the current velocity residual.
     * @return ||u^{k+1} - u^k|| / max(||u^k||, 1e-12)
     */
    [[nodiscard]] double velocityResidual() const;

    /**
     * @brief Returns the current continuity (divergence) residual.
     * @return ||∇·u^{k+1}||_2
     */
    [[nodiscard]] double continuityResidual() const;

    /**
     * @brief Returns the pressure-solver residual from the last step.
     * @return ||A p' - b||_2 from the last pressure-correction solve.
     */
    [[nodiscard]] double pressureResidual() const;

    /**
     * @brief Read-only access to the pressure field.
     * @throws std::logic_error if called before initialize().
     */
    [[nodiscard]] const Field<double>& pressure() const;

    /**
     * @brief Read-only access to the velocity field.
     * @throws std::logic_error if called before initialize().
     */
    [[nodiscard]] const Field<Eigen::Vector2d>& velocity() const;

private:
    // ── Config parameters (set in constructor, immutable after) ──────────────
    int         m_cfgNx{16};
    int         m_cfgNy{16};
    double      m_cfgLx{1.0};
    double      m_cfgLy{1.0};
    double      m_dt{0.01};
    double      m_rho{1.0};
    double      m_nu{0.01};
    double      m_tolerance{1e-6};
    double      m_alphaU{0.7};      ///< Velocity under-relaxation factor
    double      m_alphaP{0.3};      ///< Pressure under-relaxation factor
    int         m_vtkInterval{100}; ///< Write VTK every N iterations
    std::string m_outputDir{"output"};

    // ── Runtime state (valid only after initialize()) ─────────────────────────
    bool   m_initialized{false};
    double m_velResidual{1.0};      ///< ||u^{k+1} - u^k|| / max(||u^k||, 1e-12)
    double m_contResidual{1.0};     ///< ||∇·u^{k+1}||_2
    double m_pressureResidual{0.0}; ///< ||A p' - b||_2 from PressureSolver

    Mesh m_mesh;  ///< Mesh is a value member; address is stable, safe to hold refs to it.

    /// Pressure field; empty until initialize().
    std::optional<Field<double>> m_pressure;

    /// Velocity field; empty until initialize().
    std::optional<Field<Eigen::Vector2d>> m_velocity;

    BoundaryCondition m_bc;

    /// Pressure solver; empty until initialize() (requires loaded mesh).
    std::optional<PressureSolver> m_pressureSolver;

    VTKWriter m_vtkWriter;  ///< VTK output writer (stub until Milestone 2).

    /**
     * @brief Guards all methods that require prior initialization.
     *
     * @param callerName Name of the calling method for the error message.
     * @throws std::logic_error if m_initialized is false.
     */
    void checkInitialized(const char* callerName) const;
};
