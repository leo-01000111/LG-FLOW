#include "solver/MomentumSolver.hpp"

#include "solver/Discretization.hpp"

#include <Eigen/Sparse>
#include <Eigen/SparseLU>

#include <cmath>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

MomentumSolver::MomentumSolver(const Mesh& mesh)
    : m_mesh(&mesh)
{
}

void MomentumSolver::setBoundaryFaces(std::vector<VelocityBoundaryFace> velocityFaces,
                                      std::vector<PressureBoundaryFace> pressureFaces)
{
    if (velocityFaces.size() != pressureFaces.size())
        throw std::invalid_argument(
            "MomentumSolver::setBoundaryFaces: velocity and pressure face lists differ in size");

    m_velocityFaces = std::move(velocityFaces);
    m_pressureFaces = std::move(pressureFaces);

    const int Nx = m_mesh->Nx();
    const int Ny = m_mesh->Ny();
    m_westValue.assign(static_cast<std::size_t>(Ny), std::nullopt);
    m_eastValue.assign(static_cast<std::size_t>(Ny), std::nullopt);
    m_southValue.assign(static_cast<std::size_t>(Nx), std::nullopt);
    m_northValue.assign(static_cast<std::size_t>(Nx), std::nullopt);

    for (const VelocityBoundaryFace& bf : m_velocityFaces)
    {
        if (bf.kind != VelocityBCKind::DIRICHLET)
            continue;
        const auto i = static_cast<std::size_t>(bf.owner / Ny);
        const auto j = static_cast<std::size_t>(bf.owner % Ny);
        if      (bf.normal.x() < -0.5) m_westValue[j]  = bf.value;
        else if (bf.normal.x() >  0.5) m_eastValue[j]  = bf.value;
        else if (bf.normal.y() < -0.5) m_southValue[i] = bf.value;
        else                           m_northValue[i] = bf.value;
    }
}

void MomentumSolver::computeExplicitSource(const Field<Eigen::Vector2d>& uOld,
                                           const Field<double>&          pressure,
                                           double                        rho,
                                           ConvectionScheme              scheme,
                                           Field<double>&                bSrcX,
                                           Field<double>&                bSrcY) const
{
    const int  N             = m_mesh->numCells();
    const int  Nx            = m_mesh->Nx();
    const int  Ny            = m_mesh->Ny();
    const bool hasBoundaries = !m_velocityFaces.empty();

    // ── Decompose u^k into scalar components ─────────────────────────────────
    Field<double> ux(*m_mesh, 0.0);
    Field<double> uy(*m_mesh, 0.0);
    for (int c = 0; c < N; ++c)
    {
        ux[c] = uOld[c].x();
        uy[c] = uOld[c].y();
    }

    // ── Pressure gradient (explicit, from p^k) ────────────────────────────────
    // Central differencing via Gauss theorem; outlet faces use the fixed face
    // pressure when boundary faces are set.
    // Reference: Ferziger & Perić (2020) eq. 5.4.
    const Field<Eigen::Vector2d> gradP =
        hasBoundaries ? Discretization::gradient(pressure, *m_mesh, m_pressureFaces)
                      : Discretization::gradient(pressure, *m_mesh);

    // ── Explicit convective source: bSrc = -(u·∇)u - (1/ρ)∇p ────────────────
    // These terms remain explicit; only the viscous term is treated implicitly.
    if (scheme == ConvectionScheme::CENTRAL)
    {
        // Gauss face-average gradient, 2nd order central differencing, with
        // Dirichlet face velocities at walls and inlets.
        // Reference: Ferziger & Perić (2020) eq. 4.22.
        const Field<Eigen::Vector2d> gradUx =
            hasBoundaries ? Discretization::gradient(ux, *m_mesh, m_velocityFaces, 0)
                          : Discretization::gradient(ux, *m_mesh);
        const Field<Eigen::Vector2d> gradUy =
            hasBoundaries ? Discretization::gradient(uy, *m_mesh, m_velocityFaces, 1)
                          : Discretization::gradient(uy, *m_mesh);
        for (int c = 0; c < N; ++c)
        {
            bSrcX[c] = -uOld[c].dot(gradUx[c]) - (1.0 / rho) * gradP[c].x();
            bSrcY[c] = -uOld[c].dot(gradUy[c]) - (1.0 / rho) * gradP[c].y();
        }
        return;
    }

    // First-order upwind (donor-cell) differencing.
    // For φ and direction x: if u_x ≥ 0, backward diff; if u_x < 0, forward diff.
    // At a Dirichlet boundary face the donor is the face value at half-cell
    // distance; at a zero-gradient face the derivative is zero.
    // Reference: Ferziger & Perić (2020) eq. 4.24 (donor-cell).
    //
    // dx = spacing between cell centres in x; derived from uniform grid layout.
    // For a uniform grid spanning [0, Lx] with Nx cells: centre(0,j).x = dx/2.
    const double dx = 2.0 * m_mesh->getCellCenter(0, 0).x();
    const double dy = 2.0 * m_mesh->getCellCenter(0, 0).y();

    // One-sided derivative (φ_downstream − φ_upstream)/h with the boundary face
    // as upstream (or downstream) point when it is Dirichlet.
    const auto fromFace = [](const std::optional<Eigen::Vector2d>& face,
                             double phiP, int axis, double halfH, bool faceIsUpstream) {
        if (!face)
            return 0.0;
        const double phiB = (*face)[axis];
        return faceIsUpstream ? (phiP - phiB) / halfH : (phiB - phiP) / halfH;
    };

    for (int c = 0; c < N; ++c)
    {
        const int    i    = c / Ny;
        const int    j    = c % Ny;
        const auto   si   = static_cast<std::size_t>(i);
        const auto   sj   = static_cast<std::size_t>(j);
        const double ux_c = uOld[c].x();
        const double uy_c = uOld[c].y();

        double dux_dx = 0.0, dux_dy = 0.0, duy_dx = 0.0, duy_dy = 0.0;

        // x-direction upwind
        if (ux_c >= 0.0)
        {
            if (i > 0)
            {
                dux_dx = (ux[c] - ux[(i - 1) * Ny + j]) / dx;
                duy_dx = (uy[c] - uy[(i - 1) * Ny + j]) / dx;
            }
            else if (hasBoundaries)
            {
                dux_dx = fromFace(m_westValue[sj], ux[c], 0, 0.5 * dx, true);
                duy_dx = fromFace(m_westValue[sj], uy[c], 1, 0.5 * dx, true);
            }
        }
        else
        {
            if (i < Nx - 1)
            {
                dux_dx = (ux[(i + 1) * Ny + j] - ux[c]) / dx;
                duy_dx = (uy[(i + 1) * Ny + j] - uy[c]) / dx;
            }
            else if (hasBoundaries)
            {
                dux_dx = fromFace(m_eastValue[sj], ux[c], 0, 0.5 * dx, false);
                duy_dx = fromFace(m_eastValue[sj], uy[c], 1, 0.5 * dx, false);
            }
        }

        // y-direction upwind
        if (uy_c >= 0.0)
        {
            if (j > 0)
            {
                dux_dy = (ux[c] - ux[i * Ny + (j - 1)]) / dy;
                duy_dy = (uy[c] - uy[i * Ny + (j - 1)]) / dy;
            }
            else if (hasBoundaries)
            {
                dux_dy = fromFace(m_southValue[si], ux[c], 0, 0.5 * dy, true);
                duy_dy = fromFace(m_southValue[si], uy[c], 1, 0.5 * dy, true);
            }
        }
        else
        {
            if (j < Ny - 1)
            {
                dux_dy = (ux[i * Ny + (j + 1)] - ux[c]) / dy;
                duy_dy = (uy[i * Ny + (j + 1)] - uy[c]) / dy;
            }
            else if (hasBoundaries)
            {
                dux_dy = fromFace(m_northValue[si], ux[c], 0, 0.5 * dy, false);
                duy_dy = fromFace(m_northValue[si], uy[c], 1, 0.5 * dy, false);
            }
        }

        bSrcX[c] = -(ux_c * dux_dx + uy_c * dux_dy) - (1.0 / rho) * gradP[c].x();
        bSrcY[c] = -(ux_c * duy_dx + uy_c * duy_dy) - (1.0 / rho) * gradP[c].y();
    }
}

void MomentumSolver::solve(Field<Eigen::Vector2d>&       uStar,
                           const Field<Eigen::Vector2d>& uOld,
                           const Field<double>&          pressure,
                           double                        dt,
                           double                        rho,
                           double                        nu,
                           ConvectionScheme              scheme,
                           double                        alphaU)
{
    if (dt <= 0.0)
        throw std::invalid_argument("MomentumSolver::solve: dt must be > 0");
    if (rho <= 0.0)
        throw std::invalid_argument("MomentumSolver::solve: rho must be > 0");
    if (nu < 0.0)
        throw std::invalid_argument("MomentumSolver::solve: nu must be >= 0");
    if (alphaU <= 0.0 || alphaU > 1.0)
        throw std::invalid_argument("MomentumSolver::solve: alphaU must be in (0, 1]");
    if (&uOld.mesh() != m_mesh || &pressure.mesh() != m_mesh)
        throw std::invalid_argument("MomentumSolver::solve: field mesh mismatch");

    const int N  = m_mesh->numCells();
    const int Ny = m_mesh->Ny();

    // ── Explicit source: −(u·∇)u − (1/ρ)∇p ───────────────────────────────────
    Field<double> bSrcX(*m_mesh, 0.0);
    Field<double> bSrcY(*m_mesh, 0.0);
    computeExplicitSource(uOld, pressure, rho, scheme, bSrcX, bSrcY);

    // ── Build implicit viscous matrix ─────────────────────────────────────────
    // For component φ ∈ {ux, uy}:
    //   (V_P/dt + ν Σ_f a_f) φ*_P − ν a_f φ*_nb = (V_P/dt) φ^k_P + V_P bSrc_P
    //
    // a_f = A_face / |x_N − x_P|  for interior faces; Dirichlet boundary faces are
    // added below, zero-gradient boundary faces carry no viscous flux.
    // The matrix is symmetric and diagonally dominant → SPD.
    // Reference: Ferziger & Perić (2020) eq. 7.17.
    using Triplet = Eigen::Triplet<double>;
    std::vector<Triplet> triplets;
    triplets.reserve(static_cast<std::size_t>(5 * N));

    std::vector<double> diag(static_cast<std::size_t>(N), 0.0);
    Eigen::VectorXd bX(N), bY(N);

    // Pre-fill with inertia contribution (V_P/dt) and explicit source.
    for (int c = 0; c < N; ++c)
    {
        const double V                             = m_mesh->getCellVolume(c);
        diag[static_cast<std::size_t>(c)]          = V / dt;
        bX[c]                                      = (V / dt) * uOld[c].x() + V * bSrcX[c];
        bY[c]                                      = (V / dt) * uOld[c].y() + V * bSrcY[c];
    }

    // Add viscous face coupling: ν * a_f off-diagonals, accumulated on diagonal.
    const int nFaces = m_mesh->numFaces();
    for (int f = 0; f < nFaces; ++f)
    {
        const auto [o, n] = m_mesh->getNeighbors(f);
        if (n < 0)
            continue;  // boundary face: zero-gradient, no viscous coupling

        const int io  = o / Ny;
        const int jo  = o % Ny;
        const int in_ = n / Ny;
        const int jn  = n % Ny;

        const Eigen::Vector2d xo   = m_mesh->getCellCenter(io,  jo);
        const Eigen::Vector2d xn   = m_mesh->getCellCenter(in_, jn);
        const double          area = m_mesh->getFaceArea(f);
        const double          dist = (xn - xo).norm();
        const double          aF   = nu * area / dist;  // ν * A_f / d_f

        diag[static_cast<std::size_t>(o)] += aF;
        diag[static_cast<std::size_t>(n)] += aF;

        triplets.emplace_back(o, n, -aF);
        triplets.emplace_back(n, o, -aF);
    }

    // Dirichlet boundary faces (walls, inlets): viscous flux between the cell
    // centre and the face at half-cell distance d_b,
    //   ν·A·(u_b − u_P)/d_b  →  a_P += ν·A/d_b,  b += ν·A/d_b · u_b.
    // Zero-gradient faces carry no viscous flux.
    // Reference: Ferziger, Perić & Street (2020) Ch. 8 (wall boundary conditions).
    for (const VelocityBoundaryFace& bf : m_velocityFaces)
    {
        if (bf.kind != VelocityBCKind::DIRICHLET)
            continue;
        const double aB = nu * bf.area / bf.distance;
        diag[static_cast<std::size_t>(bf.owner)] += aB;
        bX[bf.owner] += aB * bf.value.x();
        bY[bf.owner] += aB * bf.value.y();
    }

    // Apply under-relaxation (Patankar 1980, eq. 6.36).
    // Replaces a_P with a_P/α_u and adds (1−α_u)/α_u · a_P · u^k to the RHS.
    // This preserves the mass conservation of u* so the pressure-corrected
    // velocity can be committed directly without a post-solve mixing step.
    if (alphaU < 1.0)
    {
        for (int c = 0; c < N; ++c)
        {
            const double aP    = diag[static_cast<std::size_t>(c)];
            const double relax = (1.0 - alphaU) / alphaU;
            bX[c] += relax * aP * uOld[c].x();
            bY[c] += relax * aP * uOld[c].y();
            diag[static_cast<std::size_t>(c)] = aP / alphaU;
        }
    }

    // Fill diagonal entries.
    for (int c = 0; c < N; ++c)
        triplets.emplace_back(c, c, diag[static_cast<std::size_t>(c)]);

    Eigen::SparseMatrix<double> A(N, N);
    A.setFromTriplets(triplets.begin(), triplets.end());
    A.makeCompressed();

    // ── Solve A x = b for ux and uy ──────────────────────────────────────────
    // SparseLU: direct solver, handles the symmetric positive-definite 5-point
    // stencil without iteration count risk.
    // Reference: Eigen documentation — SparseLU.
    //
    // Note: the factorization is shared between ux and uy because both components
    // use the same coefficient matrix A.
    Eigen::SparseLU<Eigen::SparseMatrix<double>> linSolver;
    linSolver.analyzePattern(A);
    linSolver.factorize(A);

    if (linSolver.info() != Eigen::Success)
        throw std::runtime_error(
            "MomentumSolver::solve: SparseLU factorization failed (singular matrix?)");

    const Eigen::VectorXd xSolX = linSolver.solve(bX);
    if (linSolver.info() != Eigen::Success)
        throw std::runtime_error(
            "MomentumSolver::solve: SparseLU solve failed for ux");

    const Eigen::VectorXd xSolY = linSolver.solve(bY);
    if (linSolver.info() != Eigen::Success)
        throw std::runtime_error(
            "MomentumSolver::solve: SparseLU solve failed for uy");

    // ── Verify finite solution and write to uStar ─────────────────────────────
    for (int c = 0; c < N; ++c)
    {
        if (!std::isfinite(xSolX[c]) || !std::isfinite(xSolY[c]))
            throw std::runtime_error(
                "MomentumSolver::solve: non-finite velocity at cell " +
                std::to_string(c));
        uStar[c].x() = xSolX[c];
        uStar[c].y() = xSolY[c];
    }
}

// ── Variable-viscosity overload ───────────────────────────────────────────────
// Identical to the scalar-nu overload except the face diffusion coefficient
// uses the arithmetic mean of nuEffField at owner and neighbour cells.
// Reference: Ferziger, Perić & Street (2020) §7.4 with variable ν_eff.

void MomentumSolver::solve(Field<Eigen::Vector2d>&       uStar,
                            const Field<Eigen::Vector2d>& uOld,
                            const Field<double>&          pressure,
                            const Field<double>&          nuEffField,
                            double                        dt,
                            double                        rho,
                            ConvectionScheme              scheme,
                            double                        alphaU)
{
    if (dt <= 0.0)
        throw std::invalid_argument("MomentumSolver::solve: dt must be > 0");
    if (rho <= 0.0)
        throw std::invalid_argument("MomentumSolver::solve: rho must be > 0");
    if (alphaU <= 0.0 || alphaU > 1.0)
        throw std::invalid_argument("MomentumSolver::solve: alphaU must be in (0, 1]");
    if (&uOld.mesh() != m_mesh || &pressure.mesh() != m_mesh)
        throw std::invalid_argument("MomentumSolver::solve: field mesh mismatch");

    const int N  = m_mesh->numCells();
    const int Ny = m_mesh->Ny();

    // ── Explicit source: −(u·∇)u − (1/ρ)∇p ───────────────────────────────────
    Field<double> bSrcX(*m_mesh, 0.0);
    Field<double> bSrcY(*m_mesh, 0.0);
    computeExplicitSource(uOld, pressure, rho, scheme, bSrcX, bSrcY);

    // ── Build implicit viscous matrix (variable ν_eff) ────────────────────────
    using Triplet = Eigen::Triplet<double>;
    std::vector<Triplet> triplets;
    triplets.reserve(static_cast<std::size_t>(5 * N));

    std::vector<double> diag(static_cast<std::size_t>(N), 0.0);
    Eigen::VectorXd bX(N), bY(N);

    for (int c = 0; c < N; ++c)
    {
        const double V                             = m_mesh->getCellVolume(c);
        diag[static_cast<std::size_t>(c)]          = V / dt;
        bX[c]                                      = (V / dt) * uOld[c].x() + V * bSrcX[c];
        bY[c]                                      = (V / dt) * uOld[c].y() + V * bSrcY[c];
    }

    // Face diffusion: arithmetic mean of nuEffField at owner and neighbour.
    const int nFaces = m_mesh->numFaces();
    for (int f = 0; f < nFaces; ++f)
    {
        const auto [o, n] = m_mesh->getNeighbors(f);
        if (n < 0)
            continue;

        const int io  = o / Ny;
        const int jo  = o % Ny;
        const int in_ = n / Ny;
        const int jn  = n % Ny;

        const Eigen::Vector2d xo   = m_mesh->getCellCenter(io,  jo);
        const Eigen::Vector2d xn   = m_mesh->getCellCenter(in_, jn);
        const double          area = m_mesh->getFaceArea(f);
        const double          dist = (xn - xo).norm();
        // Face-averaged effective viscosity replaces the scalar nu.
        const double          nuF  = 0.5 * (nuEffField[o] + nuEffField[n]);
        const double          aF   = nuF * area / dist;

        diag[static_cast<std::size_t>(o)] += aF;
        diag[static_cast<std::size_t>(n)] += aF;

        triplets.emplace_back(o, n, -aF);
        triplets.emplace_back(n, o, -aF);
    }

    // Dirichlet boundary faces: same half-cell wall flux as the scalar overload,
    // with ν_eff taken from the owner cell (the face has no second cell).
    for (const VelocityBoundaryFace& bf : m_velocityFaces)
    {
        if (bf.kind != VelocityBCKind::DIRICHLET)
            continue;
        const double aB = nuEffField[bf.owner] * bf.area / bf.distance;
        diag[static_cast<std::size_t>(bf.owner)] += aB;
        bX[bf.owner] += aB * bf.value.x();
        bY[bf.owner] += aB * bf.value.y();
    }

    if (alphaU < 1.0)
    {
        for (int c = 0; c < N; ++c)
        {
            const double aP    = diag[static_cast<std::size_t>(c)];
            const double relax = (1.0 - alphaU) / alphaU;
            bX[c] += relax * aP * uOld[c].x();
            bY[c] += relax * aP * uOld[c].y();
            diag[static_cast<std::size_t>(c)] = aP / alphaU;
        }
    }

    for (int c = 0; c < N; ++c)
        triplets.emplace_back(c, c, diag[static_cast<std::size_t>(c)]);

    Eigen::SparseMatrix<double> A(N, N);
    A.setFromTriplets(triplets.begin(), triplets.end());
    A.makeCompressed();

    Eigen::SparseLU<Eigen::SparseMatrix<double>> linSolver;
    linSolver.analyzePattern(A);
    linSolver.factorize(A);

    if (linSolver.info() != Eigen::Success)
        throw std::runtime_error(
            "MomentumSolver::solve (variable nu): SparseLU factorization failed");

    const Eigen::VectorXd xSolX = linSolver.solve(bX);
    if (linSolver.info() != Eigen::Success)
        throw std::runtime_error(
            "MomentumSolver::solve (variable nu): SparseLU solve failed for ux");

    const Eigen::VectorXd xSolY = linSolver.solve(bY);
    if (linSolver.info() != Eigen::Success)
        throw std::runtime_error(
            "MomentumSolver::solve (variable nu): SparseLU solve failed for uy");

    for (int c = 0; c < N; ++c)
    {
        if (!std::isfinite(xSolX[c]) || !std::isfinite(xSolY[c]))
            throw std::runtime_error(
                "MomentumSolver::solve (variable nu): non-finite velocity at cell " +
                std::to_string(c));
        uStar[c].x() = xSolX[c];
        uStar[c].y() = xSolY[c];
    }
}
