#pragma once

#include "io/VTKWriter.hpp"
#include "mesh/UnstructuredBoundaryCondition.hpp"
#include "mesh/UnstructuredField.hpp"
#include "mesh/UnstructuredMesh.hpp"
#include "solver/UnstructuredDiscretization.hpp"
#include "utils/Config.hpp"

#include <Eigen/Dense>
#include <optional>
#include <string>

/**
 * @brief SIMPLE solver for 2D incompressible Navier-Stokes on unstructured meshes.
 *
 * Implements the SIMPLE (Semi-Implicit Method for Pressure-Linked Equations)
 * algorithm adapted for unstructured FVM meshes. Mesh is loaded from a Gmsh
 * .msh v2 file specified in the config.
 *
 * SIMPLE loop order (must not be reordered without explanation):
 *   1. Momentum predictor (implicit viscous, explicit convection + pressure gradient)
 *   2. Pressure correction (Poisson equation solved with Eigen::SparseLU)
 *   3. Velocity correction (u -= (dt/ρ) ∇p')
 *   4. Apply boundary conditions
 *
 * Config keys:
 *   mesh.file          — path to Gmsh .msh file
 *   solver.dt          — time step [s] (default 0.01)
 *   solver.rho         — density [kg/m³] (default 1.0)
 *   solver.nu          — kinematic viscosity [m²/s] (default 0.01)
 *   solver.tolerance   — convergence tolerance (default 1e-6)
 *   solver.alpha_u     — velocity under-relaxation (default 0.7)
 *   solver.alpha_p     — pressure under-relaxation (default 0.3)
 *   output.vtk_interval — VTK output every N iters (default 100)
 *   output.dir         — output directory (default "output")
 *   bc.<patch>.type    — BC type (INLET, OUTLET, WALL, SYMMETRY)
 *   bc.<patch>.value_x — x-component of prescribed value
 *   bc.<patch>.value_y — y-component of prescribed value
 *
 * Reference: Patankar (1980); Ferziger, Perić & Street (2020) Ch. 7.
 */
class UnstructuredNavierStokesSolver
{
public:
    /**
     * @brief Constructs the solver and reads configuration.
     * @param config Parsed configuration object.
     * @throws std::invalid_argument on invalid parameter values.
     */
    explicit UnstructuredNavierStokesSolver(const Config& config);

    /**
     * @brief Loads mesh from Gmsh .msh file and allocates fields.
     *
     * Must be called before step(), run(), pressure(), or velocity().
     */
    void initialize();

    /**
     * @brief Performs one SIMPLE iteration.
     * @param dt Time step [s].
     * @throws std::logic_error if called before initialize().
     */
    void step(double dt);

    /**
     * @brief Runs the solver until convergence or maxIter is reached.
     * @param maxIter Maximum number of SIMPLE iterations.
     * @throws std::logic_error if called before initialize().
     */
    void run(int maxIter);

    /** @brief Returns the current velocity residual. */
    [[nodiscard]] double velocityResidual() const;

    /** @brief Returns the current continuity (divergence) residual. */
    [[nodiscard]] double continuityResidual() const;

    /**
     * @brief Read-only access to the pressure field.
     * @throws std::logic_error if called before initialize().
     */
    [[nodiscard]] const UnstructuredField<double>& pressure() const;

    /**
     * @brief Read-only access to the velocity field.
     * @throws std::logic_error if called before initialize().
     */
    [[nodiscard]] const UnstructuredField<Eigen::Vector2d>& velocity() const;

private:
    // Config parameters
    std::string m_meshFile;
    double      m_dt{0.01};
    double      m_rho{1.0};
    double      m_nu{0.01};
    double      m_tolerance{1e-6};
    double      m_alphaU{0.7};
    double      m_alphaP{0.3};
    int         m_nonorthCorrectors{2};
    int         m_vtkInterval{100};
    std::string m_outputDir{"output"};

    // Runtime state
    bool   m_initialized{false};
    double m_velResidual{1.0};
    double m_contResidual{1.0};

    UnstructuredMesh                              m_mesh;
    std::optional<UnstructuredField<double>>              m_pressure;
    std::optional<UnstructuredField<Eigen::Vector2d>>     m_velocity;
    UnstructuredBoundaryCondition                 m_bc;
    VTKWriter                                     m_vtkWriter;

    /** @brief Guards all methods that require prior initialization. */
    void checkInitialized(const char* caller) const;

    /**
     * @brief Momentum predictor: builds and solves a sparse viscous system.
     *
     * For each component φ ∈ {ux, uy}:
     *   (V/dt + ν Σ a_f) φ*_P - ν Σ a_f φ*_nb = b_P
     * where b_P includes inertia, explicit upwind convection, and pressure gradient.
     * Under-relaxation is embedded: a_P_mod = a_P/α_u, b_mod += (1-α_u)/α_u * a_P * φ^k.
     * Solved with Eigen::SparseLU (same matrix, two RHS vectors for ux/uy).
     *
     * Reference: Ferziger, Perić & Street (2020) eq. 7.17.
     *
     * @param uStar  Output predicted velocity (written in place).
     * @param uOld   Previous iteration velocity.
     * @param press  Current pressure field.
     * @param dt     Time step [s].
     */
    void solveMomentum(UnstructuredField<Eigen::Vector2d>&       uStar,
                       const UnstructuredField<Eigen::Vector2d>& uOld,
                       const UnstructuredField<double>&          press,
                       double                                    dt);

    /**
     * @brief Pressure correction: assembles Poisson system from div_RC(uStar).
     *
     * p' correction pinned to zero at OUTLET cells (or cell 0 if none).
     * After solve: p += alphaP * p'; u -= (dt/rho) * grad(p').
     * Solved with Eigen::SparseLU.
     *
     * Reference: Patankar (1980) eqs. 6.28–6.31.
     *
     * @param uStar  Predicted velocity (modified in place: velocity correction).
     * @param press  Pressure field (modified in place: pressure correction).
     * @param dt     Time step [s].
     */
    void solvePressure(UnstructuredField<Eigen::Vector2d>& uStar,
                       UnstructuredField<double>&          press,
                       double                             dt);
};
