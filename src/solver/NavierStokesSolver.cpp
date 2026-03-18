#include "solver/NavierStokesSolver.hpp"

#include "solver/Discretization.hpp"
#include "utils/Logger.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>

// ── Internal helpers ──────────────────────────────────────────────────────────

namespace
{

/// Converts a string to upper case in place.
void toUpper(std::string& s)
{
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
}

/// Parses a BC type string (case-insensitive) into a BoundaryType enum.
/// @throws std::invalid_argument on unknown type string.
BoundaryType parseBCType(const std::string& s)
{
    std::string upper = s;
    toUpper(upper);

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
    m_dt          = config.get<double>("solver.dt",            0.01);
    m_rho         = config.get<double>("solver.rho",           1.0);
    m_nu          = config.get<double>("solver.nu",            0.01);
    m_tolerance   = config.get<double>("solver.tolerance",     1e-6);
    m_alphaU      = config.get<double>("solver.alpha_u",       0.7);
    m_alphaP      = config.get<double>("solver.alpha_p",       0.3);
    m_vtkInterval = config.get<int>   ("output.vtk_interval",  100);
    m_outputDir   = config.get<std::string>("output.dir",      std::string("output"));

    // Convection scheme (default: upwind for robustness with explicit Euler).
    {
        std::string scheme = config.get<std::string>("solver.convection_scheme",
                                                     std::string("upwind"));
        toUpper(scheme);
        if (scheme == "UPWIND")
            m_convScheme = ConvectionScheme::UPWIND;
        else if (scheme == "CENTRAL")
            m_convScheme = ConvectionScheme::CENTRAL;
        else
            throw std::invalid_argument(
                "NavierStokesSolver: solver.convection_scheme must be 'upwind' or 'central'"
                " (got '" + config.get<std::string>("solver.convection_scheme", "") + "')");
    }

    // CFL safety parameters.
    m_maxCflConv = config.get<double>("solver.max_cfl_conv", 0.5);
    m_maxCflDiff = config.get<double>("solver.max_cfl_diff", 0.5);

    // Validate solver parameters — fail fast at construction rather than
    // partway through a long run.
    if (m_dt <= 0.0)
        throw std::invalid_argument(
            "NavierStokesSolver: solver.dt must be > 0 (got "
            + std::to_string(m_dt) + ")");
    if (m_rho <= 0.0)
        throw std::invalid_argument(
            "NavierStokesSolver: solver.rho must be > 0 (got "
            + std::to_string(m_rho) + ")");
    if (m_nu < 0.0)
        throw std::invalid_argument(
            "NavierStokesSolver: solver.nu must be >= 0 (got "
            + std::to_string(m_nu) + ")");
    if (m_tolerance <= 0.0)
        throw std::invalid_argument(
            "NavierStokesSolver: solver.tolerance must be > 0 (got "
            + std::to_string(m_tolerance) + ")");
    if (m_alphaU <= 0.0 || m_alphaU > 1.0)
        throw std::invalid_argument(
            "NavierStokesSolver: solver.alpha_u must be in (0, 1] (got "
            + std::to_string(m_alphaU) + ")");
    if (m_alphaP <= 0.0 || m_alphaP > 1.0)
        throw std::invalid_argument(
            "NavierStokesSolver: solver.alpha_p must be in (0, 1] (got "
            + std::to_string(m_alphaP) + ")");
    if (m_vtkInterval <= 0)
        throw std::invalid_argument(
            "NavierStokesSolver: output.vtk_interval must be > 0 (got "
            + std::to_string(m_vtkInterval) + ")");
    if (m_maxCflConv <= 0.0)
        throw std::invalid_argument(
            "NavierStokesSolver: solver.max_cfl_conv must be > 0 (got "
            + std::to_string(m_maxCflConv) + ")");
    if (m_maxCflDiff <= 0.0)
        throw std::invalid_argument(
            "NavierStokesSolver: solver.max_cfl_diff must be > 0 (got "
            + std::to_string(m_maxCflDiff) + ")");

    // Initialise m_dtEff to the configured dt; updated by each step() call.
    m_dtEff = m_dt;

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

    m_initialized     = true;
    m_velResidual     = 1.0;
    m_contResidual    = 1.0;
    m_pressureResidual = 0.0;

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

void NavierStokesSolver::step(double dt)
{
    checkInitialized("step");

    if (dt <= 0.0)
        throw std::invalid_argument(
            "NavierStokesSolver::step: dt must be > 0 (got "
            + std::to_string(dt) + ")");

    // Save u^k before any modification — needed for residual and under-relaxation.
    const Field<Eigen::Vector2d> uOld = m_velocity.value();
    const int Nc  = m_mesh.numCells();
    const int Nx  = m_cfgNx;
    const int Ny  = m_cfgNy;

    // ── CFL safety clamp ──────────────────────────────────────────────────────
    // Limit effective dt to satisfy convective and diffusive stability criteria.
    // Convective:  dt_conv = max_cfl_conv * h / ||u||_max
    // Diffusive:   dt_diff = max_cfl_diff * h² / ν
    // dt_eff = min(dt, dt_conv, dt_diff)
    // Reference: Courant-Friedrichs-Lewy (CFL) condition; von Neumann criterion.
    {
        const double dx = m_cfgLx / static_cast<double>(Nx);
        const double dy = m_cfgLy / static_cast<double>(Ny);
        const double h  = std::min(dx, dy);

        double maxVelNorm = 0.0;
        for (int c = 0; c < Nc; ++c)
            maxVelNorm = std::max(maxVelNorm, uOld[c].norm());

        const double dtConv = (maxVelNorm > 1e-10)
                                  ? m_maxCflConv * h / maxVelNorm
                                  : dt;
        const double dtDiff = (m_nu > 1e-15)
                                  ? m_maxCflDiff * h * h / m_nu
                                  : dt;
        m_dtEff = std::min({dt, dtConv, dtDiff});
    }

    // ── 1. Momentum predictor ─────────────────────────────────────────────────
    // Decompose velocity into scalar components for gradient/laplacian operators.
    Field<double> ux(m_mesh, 0.0);
    Field<double> uy(m_mesh, 0.0);
    for (int c = 0; c < Nc; ++c)
    {
        ux[c] = uOld[c].x();
        uy[c] = uOld[c].y();
    }

    const Field<double>          lapUx = Discretization::laplacian(ux, m_mesh);
    const Field<double>          lapUy = Discretization::laplacian(uy, m_mesh);
    const Field<Eigen::Vector2d> gradP = Discretization::gradient(m_pressure.value(), m_mesh);

    // u* = u^k + dt_eff * [ -(u^k·∇)u - (1/ρ) ∇p^k + ν ∇²u^k ]
    // Reference: Ferziger & Peric (2020) eq. 7.20 (explicit predictor form).
    Field<Eigen::Vector2d> uStar(m_mesh, Eigen::Vector2d::Zero());

    if (m_convScheme == ConvectionScheme::CENTRAL)
    {
        // Central differencing (Gauss face-average gradient).
        // Reference: Ferziger & Peric (2020) eq. 4.22.
        const Field<Eigen::Vector2d> gradUx = Discretization::gradient(ux, m_mesh);
        const Field<Eigen::Vector2d> gradUy = Discretization::gradient(uy, m_mesh);
        for (int c = 0; c < Nc; ++c)
        {
            const double convX = uOld[c].dot(gradUx[c]);  // (u·∇)ux — central
            const double convY = uOld[c].dot(gradUy[c]);  // (u·∇)uy — central
            uStar[c].x() = uOld[c].x() + m_dtEff * (-convX - (1.0 / m_rho) * gradP[c].x() + m_nu * lapUx[c]);
            uStar[c].y() = uOld[c].y() + m_dtEff * (-convY - (1.0 / m_rho) * gradP[c].y() + m_nu * lapUy[c]);
        }
    }
    else
    {
        // First-order upwind (donor-cell) differencing.
        // For each component φ and direction x: if u_x ≥ 0, use backward difference
        // (φ[i,j] − φ[i-1,j])/dx; if u_x < 0, forward difference.
        // At domain boundaries (no neighbour): zero-gradient closure (difference = 0).
        // Reference: Ferziger & Peric (2020) eq. 4.24 (donor-cell).
        const double dx = m_cfgLx / static_cast<double>(Nx);
        const double dy = m_cfgLy / static_cast<double>(Ny);

        for (int c = 0; c < Nc; ++c)
        {
            const int    i    = c / Ny;
            const int    j    = c % Ny;
            const double ux_c = uOld[c].x();
            const double uy_c = uOld[c].y();

            double dux_dx, dux_dy, duy_dx, duy_dy;

            // x-direction upwind
            if (ux_c >= 0.0)
            {
                dux_dx = (i > 0)      ? (ux[c] - ux[(i - 1) * Ny + j]) / dx : 0.0;
                duy_dx = (i > 0)      ? (uy[c] - uy[(i - 1) * Ny + j]) / dx : 0.0;
            }
            else
            {
                dux_dx = (i < Nx - 1) ? (ux[(i + 1) * Ny + j] - ux[c]) / dx : 0.0;
                duy_dx = (i < Nx - 1) ? (uy[(i + 1) * Ny + j] - uy[c]) / dx : 0.0;
            }

            // y-direction upwind
            if (uy_c >= 0.0)
            {
                dux_dy = (j > 0)      ? (ux[c] - ux[i * Ny + (j - 1)]) / dy : 0.0;
                duy_dy = (j > 0)      ? (uy[c] - uy[i * Ny + (j - 1)]) / dy : 0.0;
            }
            else
            {
                dux_dy = (j < Ny - 1) ? (ux[i * Ny + (j + 1)] - ux[c]) / dy : 0.0;
                duy_dy = (j < Ny - 1) ? (uy[i * Ny + (j + 1)] - uy[c]) / dy : 0.0;
            }

            const double convX = ux_c * dux_dx + uy_c * dux_dy;
            const double convY = ux_c * duy_dx + uy_c * duy_dy;

            uStar[c].x() = uOld[c].x() + m_dtEff * (-convX - (1.0 / m_rho) * gradP[c].x() + m_nu * lapUx[c]);
            uStar[c].y() = uOld[c].y() + m_dtEff * (-convY - (1.0 / m_rho) * gradP[c].y() + m_nu * lapUy[c]);
        }
    }

    // ── 2. Pressure correction ────────────────────────────────────────────────
    // PressureSolver::solve updates uStar and m_pressure in place.
    // Reference: Patankar (1980) eqs. 6.28–6.31.
    m_pressureResidual = m_pressureSolver->solve(
        uStar, m_pressure.value(), m_dtEff, m_rho, m_alphaP);

    // ── 3. Velocity under-relaxation ──────────────────────────────────────────
    // u^{k+1} = α_u · u* + (1 - α_u) · u^k
    Field<Eigen::Vector2d> uNew = uStar * m_alphaU + uOld * (1.0 - m_alphaU);

    // ── 4. Apply boundary conditions ──────────────────────────────────────────
    // BCs applied after each field update, not before (CLAUDE.md convention).
    m_bc.applyVelocity(uNew, m_mesh);
    m_bc.applyPressure(m_pressure.value(), m_mesh);

    // ── Residual tracking ─────────────────────────────────────────────────────
    // velocity residual = ||u^{k+1} - u^k|| / max(||u^k||, 1e-12)
    const Field<Eigen::Vector2d> diff = uNew + uOld * (-1.0);
    m_velResidual  = diff.norm() / std::max(uOld.norm(), 1e-12);

    // continuity residual = ||∇·u^{k+1}||_2
    m_contResidual = Discretization::divergence(uNew, m_mesh).norm();

    // Detect non-finite residuals immediately; never let them propagate silently.
    if (std::isnan(m_velResidual) || std::isinf(m_velResidual))
        throw std::runtime_error(
            "NavierStokesSolver::step: velocity residual is non-finite (NaN or Inf)");
    if (std::isnan(m_contResidual) || std::isinf(m_contResidual))
        throw std::runtime_error(
            "NavierStokesSolver::step: continuity residual is non-finite (NaN or Inf)");
    if (std::isnan(m_pressureResidual) || std::isinf(m_pressureResidual))
        throw std::runtime_error(
            "NavierStokesSolver::step: pressure residual is non-finite (NaN or Inf)");

    // Commit the updated velocity field.
    m_velocity.value() = uNew;
}

void NavierStokesSolver::run(int maxIter)
{
    checkInitialized("run");

    if (maxIter < 0)
        throw std::invalid_argument(
            "NavierStokesSolver::run: maxIter must be >= 0 (got "
            + std::to_string(maxIter) + ")");

    // Ensure output directory exists before opening files.
    namespace fs = std::filesystem;
    fs::create_directories(m_outputDir);

    const std::string histPath = m_outputDir + "/history.csv";
    std::ofstream histFile(histPath);
    if (!histFile.is_open())
        throw std::runtime_error(
            "NavierStokesSolver::run: cannot open history file: " + histPath);

    histFile << "iter,vel_residual,cont_residual,pressure_residual\n";
    histFile << std::scientific << std::setprecision(8);

    Logger::get().info(
        "NavierStokesSolver::run() - maxIter=" + std::to_string(maxIter)
        + " outputDir=" + m_outputDir);

    for (int iter = 0; iter < maxIter; ++iter)
    {
        step(m_dt);

        // Append one row to history CSV.
        histFile << (iter + 1) << ","
                 << m_velResidual  << ","
                 << m_contResidual << ","
                 << m_pressureResidual << "\n";

        // Periodic progress logging every 50 iterations.
        // Always includes dt_eff so CFL clamping is visible in the log.
        if ((iter + 1) % 50 == 0)
        {
            std::string msg =
                "iter=" + std::to_string(iter + 1)
                + " vel="  + std::to_string(m_velResidual)
                + " cont=" + std::to_string(m_contResidual)
                + " p="    + std::to_string(m_pressureResidual)
                + " dt_eff=" + std::to_string(m_dtEff);
            if (m_dtEff < m_dt * (1.0 - 1e-6))
                msg += " [CFL-clamped]";
            Logger::get().info(msg);
        }

        // VTK snapshot every vtk_interval iterations.
        if ((iter + 1) % m_vtkInterval == 0)
        {
            const std::string vtkFile =
                m_outputDir + "/iter_" + std::to_string(iter + 1) + ".vtu";
            m_vtkWriter.write(vtkFile, m_mesh,
                              m_pressure.value(), m_velocity.value());
        }

        // Convergence: both velocity residual AND continuity residual must
        // fall below tolerance (checking one alone is insufficient).
        if (m_velResidual < m_tolerance && m_contResidual < m_tolerance)
        {
            Logger::get().info(
                "Converged at iteration " + std::to_string(iter + 1)
                + " (vel=" + std::to_string(m_velResidual)
                + " cont=" + std::to_string(m_contResidual) + ")");
            return;
        }
    }

    if (maxIter > 0)
        Logger::get().warn(
            "run() reached maxIter=" + std::to_string(maxIter)
            + " without convergence");
}

double NavierStokesSolver::residual() const
{
    return m_velResidual;
}

double NavierStokesSolver::velocityResidual() const
{
    return m_velResidual;
}

double NavierStokesSolver::continuityResidual() const
{
    return m_contResidual;
}

double NavierStokesSolver::pressureResidual() const
{
    return m_pressureResidual;
}

double NavierStokesSolver::dtEffective() const
{
    return m_dtEff;
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
