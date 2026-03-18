#include "solver/NavierStokesSolver.hpp"

#include "utils/Logger.hpp"

NavierStokesSolver::NavierStokesSolver(const Config& config)
    : m_mesh()
    , m_pressure(m_mesh)    // bound to m_mesh; re-initialised in initialize()
    , m_velocity(m_mesh)
    , m_bc()
    , m_pressureSolver(m_mesh)
{
    // STUB: Read parameters from config.
    // Will be fully implemented in Milestone 3.
    m_dt        = config.get<double>("solver.dt",        0.01);
    m_rho       = config.get<double>("solver.rho",       1.0);
    m_nu        = config.get<double>("solver.nu",        0.01);
    m_tolerance = config.get<double>("solver.tolerance", 1e-6);
}

void NavierStokesSolver::initialize()
{
    // STUB: Build mesh, allocate fields, set up BCs, zero-initialise.
    // Will be implemented in Milestone 2 (mesh) and Milestone 3 (BCs).
    Logger::get().info("NavierStokesSolver::initialize() — stub");
}

void NavierStokesSolver::step(double /*dt*/)
{
    // STUB: SIMPLE iteration.
    // Order (must not be reordered without explanation):
    //   1. Momentum predictor
    //   2. m_pressureSolver.solve(...)
    //   3. Velocity correction
    //   4. m_bc.applyVelocity / applyPressure
    // Will be implemented in Milestone 3.
}

void NavierStokesSolver::run(int maxIter)
{
    Logger::get().info("NavierStokesSolver::run() — stub, maxIter = "
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
    return m_pressure;
}

const Field<Eigen::Vector2d>& NavierStokesSolver::velocity() const
{
    return m_velocity;
}
