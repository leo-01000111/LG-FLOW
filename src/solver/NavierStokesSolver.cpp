#include "solver/NavierStokesSolver.hpp"

#include "utils/Logger.hpp"

#include <stdexcept>
#include <string>

NavierStokesSolver::NavierStokesSolver(const Config& config)
{
    // Read all parameters from config; use safe defaults so that a
    // default-constructed (empty) Config is still valid for testing.
    m_cfgNx     = config.get<int>   ("mesh.Nx",          16);
    m_cfgNy     = config.get<int>   ("mesh.Ny",          16);
    m_cfgLx     = config.get<double>("mesh.Lx",          1.0);
    m_cfgLy     = config.get<double>("mesh.Ly",          1.0);
    m_dt        = config.get<double>("solver.dt",         0.01);
    m_rho       = config.get<double>("solver.rho",        1.0);
    m_nu        = config.get<double>("solver.nu",         0.01);
    m_tolerance = config.get<double>("solver.tolerance",  1e-6);
}

void NavierStokesSolver::initialize()
{
    // Build the mesh from config dimensions. Throws std::invalid_argument on
    // bad dimensions (zero, negative) propagated from Mesh::load().
    m_mesh.load(m_cfgNx, m_cfgNy, m_cfgLx, m_cfgLy);

    // Allocate fields bound to the now-loaded mesh, zero-initialised.
    m_pressure.emplace(m_mesh, 0.0);
    m_velocity.emplace(m_mesh, Eigen::Vector2d::Zero());

    // Construct pressure solver bound to the loaded mesh.
    m_pressureSolver.emplace(m_mesh);

    m_initialized = true;
    m_residual    = 1.0;  // reset residual history

    Logger::get().info(
        "NavierStokesSolver::initialize() - stub, mesh "
        + std::to_string(m_cfgNx) + "x" + std::to_string(m_cfgNy));
}

void NavierStokesSolver::checkInitialized(const char* callerName) const
{
    if (!m_initialized)
        throw std::logic_error(
            std::string("NavierStokesSolver::") + callerName +
            "() called before initialize() -- call initialize() first");
}

void NavierStokesSolver::step(double /*dt*/)
{
    checkInitialized("step");

    // STUB: SIMPLE iteration.
    // Order (must not be reordered without explanation):
    //   1. Momentum predictor
    //   2. m_pressureSolver->solve(...)
    //   3. Velocity correction
    //   4. m_bc.applyVelocity / applyPressure
    // Will be implemented in Milestone 3.
}

void NavierStokesSolver::run(int maxIter)
{
    checkInitialized("run");

    Logger::get().info(
        "NavierStokesSolver::run() - stub, maxIter = "
        + std::to_string(maxIter));

    for (int iter = 0; iter < maxIter; ++iter)
    {
        step(m_dt);

        if (m_residual < m_tolerance)
        {
            Logger::get().info("Converged at iteration " + std::to_string(iter));
            return;
        }
    }

    Logger::get().warn("run() reached maxIter without convergence (stub)");
}

double NavierStokesSolver::residual() const
{
    return m_residual;
}

const Field<double>& NavierStokesSolver::pressure() const
{
    checkInitialized("pressure");
    return m_pressure.value();
}

const Field<Eigen::Vector2d>& NavierStokesSolver::velocity() const
{
    checkInitialized("velocity");
    return m_velocity.value();
}
