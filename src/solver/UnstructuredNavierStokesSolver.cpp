#include "solver/UnstructuredNavierStokesSolver.hpp"

#include "io/GmshReader.hpp"
#include "utils/Logger.hpp"

#include <Eigen/Sparse>
#include <Eigen/SparseLU>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <vector>

// ── Local helpers ─────────────────────────────────────────────────────────────

namespace
{

void toUpperInPlace(std::string& s)
{
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
}

BoundaryType parseBCType(const std::string& s)
{
    std::string upper = s;
    toUpperInPlace(upper);

    if (upper == "INLET")           return BoundaryType::INLET;
    if (upper == "PARABOLIC_INLET") return BoundaryType::PARABOLIC_INLET;
    if (upper == "OUTLET")          return BoundaryType::OUTLET;
    if (upper == "WALL")            return BoundaryType::WALL;
    if (upper == "SYMMETRY")        return BoundaryType::SYMMETRY;

    throw std::invalid_argument(
        "UnstructuredNavierStokesSolver: unknown BC type '" + s + "'");
}

} // namespace

// ── Constructor ───────────────────────────────────────────────────────────────

UnstructuredNavierStokesSolver::UnstructuredNavierStokesSolver(const Config& config)
{
    m_meshFile   = config.get<std::string>("mesh.file", std::string("mesh.msh"));
    m_dt         = config.get<double>("solver.dt",         0.01);
    m_rho        = config.get<double>("solver.rho",        1.0);
    m_nu         = config.get<double>("solver.nu",         0.01);
    m_tolerance  = config.get<double>("solver.tolerance",  1e-6);
    m_alphaU     = config.get<double>("solver.alpha_u",    0.7);
    m_alphaP     = config.get<double>("solver.alpha_p",    0.3);
    m_nonorthCorrectors = config.get<int>("solver.nonorth_correctors", 2);
    m_vtkInterval = config.get<int>("output.vtk_interval", 100);
    m_outputDir   = config.get<std::string>("output.dir",  std::string("output"));

    if (m_dt <= 0.0)
        throw std::invalid_argument(
            "UnstructuredNavierStokesSolver: solver.dt must be > 0");
    if (m_rho <= 0.0)
        throw std::invalid_argument(
            "UnstructuredNavierStokesSolver: solver.rho must be > 0");
    if (m_nu < 0.0)
        throw std::invalid_argument(
            "UnstructuredNavierStokesSolver: solver.nu must be >= 0");
    if (m_alphaU <= 0.0 || m_alphaU > 1.0)
        throw std::invalid_argument(
            "UnstructuredNavierStokesSolver: solver.alpha_u must be in (0, 1]");
    if (m_alphaP <= 0.0 || m_alphaP > 1.0)
        throw std::invalid_argument(
            "UnstructuredNavierStokesSolver: solver.alpha_p must be in (0, 1]");

    // Parse boundary conditions for all patches listed in config.
    // Keys: bc.<patchname>.type, bc.<patchname>.value_x, bc.<patchname>.value_y
    // We scan common names; actual patch validation happens at run time vs mesh patches.
    //
    // Strategy: iterate all config keys to find bc.*.type entries.
    // Config does not expose a key-enumeration API, so we probe a fixed set of
    // common patch names plus any that may be in the mesh.  Unrecognised bc keys
    // are silently skipped.
    //
    // For channel_unstruct: inlet, outlet, wall_bottom, wall_top
    // This list covers typical use cases; extend as needed.
    const std::vector<std::string> candidatePatches{
        "inlet", "outlet", "wall", "wall_bottom", "wall_top",
        "left",  "right",  "bottom", "top", "symmetry",
        "inlet_left", "outlet_right", "farfield", "moving_wall"
    };

    for (const auto& patchName : candidatePatches)
    {
        const std::string typeKey = "bc." + patchName + ".type";
        if (!config.has(typeKey))
            continue;

        const BoundaryType bt = parseBCType(config.get<std::string>(typeKey));

        double vx = 0.0;
        if (config.has("bc." + patchName + ".value_x"))
            vx = config.get<double>("bc." + patchName + ".value_x");
        else if (config.has("bc." + patchName + ".value"))
            vx = config.get<double>("bc." + patchName + ".value");

        double vy = 0.0;
        if (config.has("bc." + patchName + ".value_y"))
            vy = config.get<double>("bc." + patchName + ".value_y");

        m_bc.addPatch(patchName, {bt, Eigen::Vector2d(vx, vy)});
    }
}

// ── initialize ────────────────────────────────────────────────────────────────

void UnstructuredNavierStokesSolver::initialize()
{
    GmshReader reader;
    reader.read(m_meshFile, m_mesh);

    m_pressure.emplace(m_mesh, 0.0);
    m_velocity.emplace(m_mesh, Eigen::Vector2d::Zero());

    m_bc.applyVelocity(m_velocity.value(), m_mesh);
    m_bc.applyPressure(m_pressure.value(), m_mesh);

    m_initialized  = true;
    m_velResidual  = 1.0;
    m_contResidual = 1.0;

    Logger::get().info(
        "UnstructuredNavierStokesSolver::initialize() - "
        + std::to_string(m_mesh.numCells()) + " cells, "
        + std::to_string(m_mesh.numFaces()) + " faces");
}

// ── checkInitialized ──────────────────────────────────────────────────────────

void UnstructuredNavierStokesSolver::checkInitialized(const char* caller) const
{
    if (!m_initialized)
        throw std::logic_error(
            std::string("UnstructuredNavierStokesSolver::") + caller +
            "() called before initialize() -- call initialize() first");
}

// ── solveMomentum ─────────────────────────────────────────────────────────────
//
// Assembles and solves the implicit viscous momentum equation.
// Diagonal = V/dt + ν * Σ_f a_f (diffusion face coefficients).
// Off-diagonal = -ν * a_f (viscous coupling between P and N).
// RHS: inertia + explicit upwind convection source - (1/ρ) * V * ∇p.
// Under-relaxation embedded (Patankar 1980, eq. 6.36).
// Reference: Ferziger, Perić & Street (2020) eq. 7.17.

void UnstructuredNavierStokesSolver::solveMomentum(
    UnstructuredField<Eigen::Vector2d>&       uStar,
    const UnstructuredField<Eigen::Vector2d>& uOld,
    const UnstructuredField<double>&          press,
    double                                    dt)
{
    const int N = m_mesh.numCells();

    // Pressure gradient (explicit).
    const UnstructuredField<Eigen::Vector2d> gradP =
        UnstructuredDiscretization::gradient(press, m_mesh);

    using Triplet = Eigen::Triplet<double>;
    std::vector<Triplet> triplets;
    triplets.reserve(static_cast<std::size_t>(5 * N));

    std::vector<double> diag(static_cast<std::size_t>(N), 0.0);
    Eigen::VectorXd bX(N), bY(N);

    // Inertia contribution: V/dt on diagonal; (V/dt)*u^k on RHS.
    for (int ci = 0; ci < N; ++ci)
    {
        const double V  = m_mesh.cell(CellId{ci}).volume;
        diag[static_cast<std::size_t>(ci)] = V / dt;
        bX[ci] = (V / dt) * uOld[CellId{ci}].x() - V * (1.0 / m_rho) * gradP[CellId{ci}].x();
        bY[ci] = (V / dt) * uOld[CellId{ci}].y() - V * (1.0 / m_rho) * gradP[CellId{ci}].y();
    }

    // Explicit upwind convection source.
    // For each cell, iterate its faces. Use upwind-weighted value based on
    // face-normal velocity sign (inflow vs outflow from owner perspective).
    for (int ci = 0; ci < N; ++ci)
    {
        const MeshCell& c     = m_mesh.cell(CellId{ci});
        const double    V     = c.volume;
        double          srcX  = 0.0;
        double          srcY  = 0.0;

        for (FaceId fid : c.faces)
        {
            const MeshFace& f = m_mesh.face(fid);

            // Outward normal sign: if this cell is the owner, normal points away;
            // if this cell is the neighbour, normal points toward this cell.
            const bool isOwner = (toInt(f.owner) == ci);
            const Eigen::Vector2d n_out = isOwner ? f.normal : -f.normal;

            const Eigen::Vector2d u_P = uOld[CellId{ci}];

            // Face velocity for upwind: determine inflow or outflow using the owner's perspective.
            // Face normal is always stored from owner's outward perspective.
            const Eigen::Vector2d u_face_avg =
                (toInt(f.neighbour) < 0)
                    ? uOld[f.owner]
                    : 0.5 * (uOld[f.owner] + uOld[f.neighbour]);

            const double dotProduct = u_face_avg.dot(f.normal);

            // Upwind pick: if flow is out of owner (dot > 0), use owner value;
            // if flow is into owner (dot < 0), use neighbour value.
            Eigen::Vector2d u_upwind;
            if (toInt(f.neighbour) < 0)
            {
                u_upwind = uOld[f.owner];
            }
            else if (dotProduct >= 0.0)
            {
                u_upwind = uOld[f.owner];
            }
            else
            {
                u_upwind = uOld[f.neighbour];
            }

            // Flux contribution to this cell: inward flux is negative for outward normal.
            const double flux_sign = isOwner ? 1.0 : -1.0;
            const double conv_flux = flux_sign * u_face_avg.dot(f.normal) * f.area;

            // Convective source: -(u·∇)u approximated as face flux * upwind value / V.
            srcX -= conv_flux * u_upwind.x() / V;
            srcY -= conv_flux * u_upwind.y() / V;

            // Suppress unused variable warning (u_P captured for context but not used directly).
            (void)u_P;
            (void)n_out;
        }

        bX[ci] += m_mesh.cell(CellId{ci}).volume * srcX;
        bY[ci] += m_mesh.cell(CellId{ci}).volume * srcY;
    }

    // Viscous diffusion (implicit): for each interior face, add ν*a_f coupling.
    // a_f = A_face / |x_N - x_P|  — face diffusion coefficient.
    // Reference: Ferziger, Perić & Street (2020) eq. 7.17.
    const int nFaces = m_mesh.numFaces();
    for (int fi = 0; fi < nFaces; ++fi)
    {
        const MeshFace& f = m_mesh.face(FaceId{fi});
        if (toInt(f.neighbour) < 0)
            continue;  // boundary face: zero-gradient closure, no off-diagonal entry

        const int iP = toInt(f.owner);
        const int iN = toInt(f.neighbour);

        const Eigen::Vector2d& xP = m_mesh.cell(CellId{iP}).centre;
        const Eigen::Vector2d& xN = m_mesh.cell(CellId{iN}).centre;
        const double dist = (xN - xP).norm();
        const double aF   = m_nu * f.area / dist;

        diag[static_cast<std::size_t>(iP)] += aF;
        diag[static_cast<std::size_t>(iN)] += aF;

        triplets.emplace_back(iP, iN, -aF);
        triplets.emplace_back(iN, iP, -aF);
    }

    // Under-relaxation (Patankar 1980, eq. 6.36).
    for (int ci = 0; ci < N; ++ci)
    {
        const double aP    = diag[static_cast<std::size_t>(ci)];
        const double relax = (1.0 - m_alphaU) / m_alphaU;
        bX[ci] += relax * aP * uOld[CellId{ci}].x();
        bY[ci] += relax * aP * uOld[CellId{ci}].y();
        diag[static_cast<std::size_t>(ci)] = aP / m_alphaU;
    }

    // Fill diagonal entries.
    for (int ci = 0; ci < N; ++ci)
        triplets.emplace_back(ci, ci, diag[static_cast<std::size_t>(ci)]);

    Eigen::SparseMatrix<double> A(N, N);
    A.setFromTriplets(triplets.begin(), triplets.end());
    A.makeCompressed();

    Eigen::SparseLU<Eigen::SparseMatrix<double>> linSolver;
    linSolver.analyzePattern(A);
    linSolver.factorize(A);

    if (linSolver.info() != Eigen::Success)
        throw std::runtime_error(
            "UnstructuredNavierStokesSolver::solveMomentum: SparseLU factorization failed");

    const Eigen::VectorXd solX = linSolver.solve(bX);
    if (linSolver.info() != Eigen::Success)
        throw std::runtime_error(
            "UnstructuredNavierStokesSolver::solveMomentum: SparseLU solve failed for ux");

    const Eigen::VectorXd solY = linSolver.solve(bY);
    if (linSolver.info() != Eigen::Success)
        throw std::runtime_error(
            "UnstructuredNavierStokesSolver::solveMomentum: SparseLU solve failed for uy");

    for (int ci = 0; ci < N; ++ci)
    {
        if (!std::isfinite(solX[ci]) || !std::isfinite(solY[ci]))
            throw std::runtime_error(
                "UnstructuredNavierStokesSolver::solveMomentum: non-finite velocity at cell "
                + std::to_string(ci));
        uStar[CellId{ci}].x() = solX[ci];
        uStar[CellId{ci}].y() = solY[ci];
    }
}

// ── solvePressure ─────────────────────────────────────────────────────────────
//
// Assembles and solves the pressure-correction Poisson equation.
// a_f = A_face / |x_N - x_P| for each interior face.
// RHS: b_c = -(ρ/dt) * div_RC(u*) * V_c
// Reference cells: OUTLET cells pinned to p' = 0 (or cell 0 if no outlet).
// Reference: Patankar (1980) eqs. 6.28–6.31.

void UnstructuredNavierStokesSolver::solvePressure(
    UnstructuredField<Eigen::Vector2d>& uStar,
    UnstructuredField<double>&          press,
    double                             dt)
{
    const int N = m_mesh.numCells();

    // Divergence of predicted velocity (Rhie-Chow corrected).
    const UnstructuredField<double> divU =
        UnstructuredDiscretization::divergenceRhieChow(uStar, press, m_mesh, dt, m_rho);

    Eigen::VectorXd b(N);
    for (int ci = 0; ci < N; ++ci)
    {
        const double V = m_mesh.cell(CellId{ci}).volume;
        b[ci] = -(m_rho / dt) * divU[CellId{ci}] * V;
    }

    // Identify reference cells (OUTLET or fallback to cell 0).
    const std::vector<CellId> refCellIds =
        m_bc.collectCellsOfType(BoundaryType::OUTLET, m_mesh);

    std::vector<bool> isRef(static_cast<std::size_t>(N), false);
    if (refCellIds.empty())
    {
        isRef[0] = true;
    }
    else
    {
        for (CellId cid : refCellIds)
        {
            const int idx = toInt(cid);
            if (idx >= 0 && idx < N)
                isRef[static_cast<std::size_t>(idx)] = true;
        }
    }

    using Triplet = Eigen::Triplet<double>;
    std::vector<Triplet> triplets;
    triplets.reserve(static_cast<std::size_t>(5 * N));

    std::vector<double> diag(static_cast<std::size_t>(N), 0.0);

    const int nFaces = m_mesh.numFaces();
    for (int fi = 0; fi < nFaces; ++fi)
    {
        const MeshFace& f = m_mesh.face(FaceId{fi});
        if (toInt(f.neighbour) < 0)
            continue;

        const int iP = toInt(f.owner);
        const int iN = toInt(f.neighbour);

        const Eigen::Vector2d& xP = m_mesh.cell(CellId{iP}).centre;
        const Eigen::Vector2d& xN = m_mesh.cell(CellId{iN}).centre;
        const double dist = (xN - xP).norm();
        const double aF   = f.area / dist;

        diag[static_cast<std::size_t>(iP)] += aF;
        diag[static_cast<std::size_t>(iN)] += aF;

        if (!isRef[static_cast<std::size_t>(iP)])
            triplets.emplace_back(iP, iN, -aF);
        if (!isRef[static_cast<std::size_t>(iN)])
            triplets.emplace_back(iN, iP, -aF);
    }

    // Apply reference cell identity rows.
    for (int ci = 0; ci < N; ++ci)
    {
        if (isRef[static_cast<std::size_t>(ci)])
        {
            triplets.emplace_back(ci, ci, 1.0);
            b[ci] = 0.0;
        }
        else
        {
            triplets.emplace_back(ci, ci, diag[static_cast<std::size_t>(ci)]);
        }
    }

    Eigen::SparseMatrix<double> A(N, N);
    A.setFromTriplets(triplets.begin(), triplets.end());
    A.makeCompressed();

    Eigen::SparseLU<Eigen::SparseMatrix<double>> directSolver;
    directSolver.analyzePattern(A);
    directSolver.factorize(A);

    if (directSolver.info() != Eigen::Success)
        throw std::runtime_error(
            "UnstructuredNavierStokesSolver::solvePressure: SparseLU factorization failed");

    const Eigen::VectorXd x = directSolver.solve(b);

    if (directSolver.info() != Eigen::Success)
        throw std::runtime_error(
            "UnstructuredNavierStokesSolver::solvePressure: SparseLU solve failed");

    for (int ci = 0; ci < N; ++ci)
    {
        if (!std::isfinite(x[ci]))
            throw std::runtime_error(
                "UnstructuredNavierStokesSolver::solvePressure: non-finite p' at cell "
                + std::to_string(ci));
    }

    // Correct pressure: p += alphaP * p'.
    for (int ci = 0; ci < N; ++ci)
        press[CellId{ci}] += m_alphaP * x[ci];

    // Correct velocity: u -= (dt/rho) * grad(p').
    UnstructuredField<double> pPrime(m_mesh, 0.0);
    for (int ci = 0; ci < N; ++ci)
        pPrime[CellId{ci}] = x[ci];

    const UnstructuredField<Eigen::Vector2d> gradPPrime =
        UnstructuredDiscretization::gradient(pPrime, m_mesh);

    const double coeff = dt / m_rho;
    for (int ci = 0; ci < N; ++ci)
        uStar[CellId{ci}] -= coeff * gradPPrime[CellId{ci}];
}

// ── step ─────────────────────────────────────────────────────────────────────

void UnstructuredNavierStokesSolver::step(double dt)
{
    checkInitialized("step");

    if (dt <= 0.0)
        throw std::invalid_argument(
            "UnstructuredNavierStokesSolver::step: dt must be > 0");

    const UnstructuredField<Eigen::Vector2d> uOld = m_velocity.value();

    // 1. Momentum predictor.
    UnstructuredField<Eigen::Vector2d> uStar(m_mesh, Eigen::Vector2d::Zero());
    solveMomentum(uStar, uOld, m_pressure.value(), dt);

    // 1b. Apply BCs to u* before pressure correction.
    m_bc.applyVelocity(uStar, m_mesh);

    // Continuity residual of u* (before pressure correction).
    m_contResidual =
        UnstructuredDiscretization::divergenceRhieChow(
            uStar, m_pressure.value(), m_mesh, dt, m_rho).norm();

    // 2. Pressure correction.
    solvePressure(uStar, m_pressure.value(), dt);

    // 3. Commit corrected velocity.
    UnstructuredField<Eigen::Vector2d> uNew = uStar;

    // 4. Apply boundary conditions.
    m_bc.applyVelocity(uNew, m_mesh);
    m_bc.applyPressure(m_pressure.value(), m_mesh);

    // Velocity residual.
    const UnstructuredField<Eigen::Vector2d> diff = uNew + uOld * (-1.0);
    m_velResidual = diff.norm() / std::max(uOld.norm(), 1e-12);

    if (std::isnan(m_velResidual) || std::isinf(m_velResidual))
        throw std::runtime_error(
            "UnstructuredNavierStokesSolver::step: velocity residual is non-finite");
    if (std::isnan(m_contResidual) || std::isinf(m_contResidual))
        throw std::runtime_error(
            "UnstructuredNavierStokesSolver::step: continuity residual is non-finite");

    m_velocity.value() = uNew;
}

// ── run ───────────────────────────────────────────────────────────────────────

void UnstructuredNavierStokesSolver::run(int maxIter)
{
    checkInitialized("run");

    if (maxIter < 0)
        throw std::invalid_argument(
            "UnstructuredNavierStokesSolver::run: maxIter must be >= 0");

    namespace fs = std::filesystem;
    fs::create_directories(m_outputDir);

    const std::string histPath = m_outputDir + "/history.csv";
    std::ofstream histFile(histPath);
    if (!histFile.is_open())
        throw std::runtime_error(
            "UnstructuredNavierStokesSolver::run: cannot open history file: " + histPath);

    histFile << "iter,vel_residual,cont_residual\n";
    histFile << std::scientific << std::setprecision(8);

    Logger::get().info(
        "UnstructuredNavierStokesSolver::run() - maxIter=" + std::to_string(maxIter));

    for (int iter = 0; iter < maxIter; ++iter)
    {
        step(m_dt);

        histFile << (iter + 1) << ","
                 << m_velResidual  << ","
                 << m_contResidual << "\n";

        if ((iter + 1) % 50 == 0)
        {
            Logger::get().info(
                "iter=" + std::to_string(iter + 1)
                + " vel="  + std::to_string(m_velResidual)
                + " cont=" + std::to_string(m_contResidual));
        }

        if ((iter + 1) % m_vtkInterval == 0)
        {
            const std::string vtkFile =
                m_outputDir + "/iter_" + std::to_string(iter + 1) + ".vtk";
            m_vtkWriter.write(vtkFile, m_mesh,
                              m_pressure.value(), m_velocity.value());
        }

        if (m_velResidual < m_tolerance)
        {
            Logger::get().info(
                "Converged at iteration " + std::to_string(iter + 1)
                + " (vel=" + std::to_string(m_velResidual) + ")");
            return;
        }
    }

    if (maxIter > 0)
        Logger::get().warn(
            "run() reached maxIter=" + std::to_string(maxIter)
            + " without convergence");
}

// ── Query methods ─────────────────────────────────────────────────────────────

double UnstructuredNavierStokesSolver::velocityResidual() const
{
    return m_velResidual;
}

double UnstructuredNavierStokesSolver::continuityResidual() const
{
    return m_contResidual;
}

const UnstructuredField<double>& UnstructuredNavierStokesSolver::pressure() const
{
    checkInitialized("pressure");
    return m_pressure.value();
}

const UnstructuredField<Eigen::Vector2d>& UnstructuredNavierStokesSolver::velocity() const
{
    checkInitialized("velocity");
    return m_velocity.value();
}
