#include "solver/NavierStokesSolver.hpp"

#include "utils/Logger.hpp"

#include <algorithm>
#include <cctype>
#include <stdexcept>
#include <string>

// ── Internal helpers ──────────────────────────────────────────────────────────

namespace
{

/// Parses a BC type string (case-insensitive) into a BoundaryType enum.
/// @throws std::invalid_argument on unknown type string.
BoundaryType parseBCType(const std::string& s)
{
    std::string upper = s;
    std::transform(upper.begin(), upper.end(), upper.begin(),
                   [](unsigned char c) { return static_cast<char>(std::toupper(c)); });

    if (upper == "INLET")    return BoundaryType::INLET;
    if (upper == "OUTLET")   return BoundaryType::OUTLET;
    if (upper == "WALL")     return BoundaryType::WALL;
    if (upper == "SYMMETRY") return BoundaryType::SYMMETRY;

    throw std::invalid_argument(
        "NavierStokesSolver: unknown BC type '" + s + "'"
        " (expected INLET, OUTLET, WALL, or SYMMETRY)");
}

/// Builds and returns a BoundaryCondition populated from config bc.* keys.
/// Sides without a bc.<side>.type key are silently skipped.
/// Value parsing: x = bc.<side>.value_x > bc.<side>.value > 0.0
///                y = bc.<side>.value_y > 0.0
BoundaryCondition buildBoundaryConditions(const Config& config)
{
    BoundaryCondition bc;

    for (const std::string& side : std::initializer_list<std::string>{"left", "right", "bottom", "top"})
    {
        const std::string typeKey = "bc." + side + ".type";
        if (!config.has(typeKey))
            continue;

        const BoundaryType bt = parseBCType(config.get<std::string>(typeKey));

        const std::string vxKey = "bc." + side + ".value_x";
        const std::string vKey  = "bc." + side + ".value";
        const std::string vyKey = "bc." + side + ".value_y";

        double vx = 0.0;
        if (config.has(vxKey))
            vx = config.get<double>(vxKey);
        else if (config.has(vKey))
            vx = config.get<double>(vKey);

        double vy = 0.0;
        if (config.has(vyKey))
            vy = config.get<double>(vyKey);

        bc.addPatch(side, {bt, Eigen::Vector2d(vx, vy)});
    }

    return bc;
}

}  // namespace

// ── NavierStokesSolver ────────────────────────────────────────────────────────

NavierStokesSolver::NavierStokesSolver(const Config& config)
{
    // Read mesh and solver parameters from config with safe defaults,
    // so that a default-constructed (empty) Config is still valid for testing.
    m_cfgNx     = config.get<int>   ("mesh.Nx",          16);
    m_cfgNy     = config.get<int>   ("mesh.Ny",          16);
    m_cfgLx     = config.get<double>("mesh.Lx",          1.0);
    m_cfgLy     = config.get<double>("mesh.Ly",          1.0);
    m_dt        = config.get<double>("solver.dt",         0.01);
    m_rho       = config.get<double>("solver.rho",        1.0);
    m_nu        = config.get<double>("solver.nu",         0.01);
    m_tolerance = config.get<double>("solver.tolerance",  1e-6);

    // Parse boundary condition config; throws std::invalid_argument on unknown type.
    m_bc = buildBoundaryConditions(config);
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

    // Apply boundary conditions to set initial field values at boundaries.
    // Boundary conditions are applied after each field update in the SIMPLE loop;
    // applying once here sets correct initial values for the first iteration.
    m_bc.applyVelocity(m_velocity.value(), m_mesh);
    m_bc.applyPressure(m_pressure.value(), m_mesh);

    m_initialized = true;
    m_residual    = 1.0;  // reset residual history

    Logger::get().info(
        "NavierStokesSolver::initialize() - mesh "
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
