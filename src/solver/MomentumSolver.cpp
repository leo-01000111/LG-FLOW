#include "solver/MomentumSolver.hpp"

#include "solver/Discretization.hpp"

#include <Eigen/Sparse>
#include <Eigen/SparseLU>

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

MomentumSolver::MomentumSolver(const Mesh& mesh)
    : m_mesh(&mesh)
{
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
    const int Nx = m_mesh->Nx();
    const int Ny = m_mesh->Ny();

    // ── Decompose u^k into scalar components ─────────────────────────────────
    Field<double> ux(*m_mesh, 0.0);
    Field<double> uy(*m_mesh, 0.0);
    for (int c = 0; c < N; ++c)
    {
        ux[c] = uOld[c].x();
        uy[c] = uOld[c].y();
    }

    // ── Pressure gradient (explicit, from p^k) ────────────────────────────────
    // Central differencing via Gauss theorem.
    // Reference: Ferziger & Perić (2020) eq. 5.4.
    const Field<Eigen::Vector2d> gradP = Discretization::gradient(pressure, *m_mesh);

    // ── Explicit convective source: bSrc = -(u·∇)u - (1/ρ)∇p ────────────────
    // These terms remain explicit; only the viscous term is treated implicitly.
    Field<double> bSrcX(*m_mesh, 0.0);
    Field<double> bSrcY(*m_mesh, 0.0);

    if (scheme == ConvectionScheme::CENTRAL)
    {
        // Gauss face-average gradient, 2nd order central differencing.
        // Reference: Ferziger & Perić (2020) eq. 4.22.
        const Field<Eigen::Vector2d> gradUx = Discretization::gradient(ux, *m_mesh);
        const Field<Eigen::Vector2d> gradUy = Discretization::gradient(uy, *m_mesh);
        for (int c = 0; c < N; ++c)
        {
            bSrcX[c] = -uOld[c].dot(gradUx[c]) - (1.0 / rho) * gradP[c].x();
            bSrcY[c] = -uOld[c].dot(gradUy[c]) - (1.0 / rho) * gradP[c].y();
        }
    }
    else
    {
        // First-order upwind (donor-cell) differencing.
        // For φ and direction x: if u_x ≥ 0, backward diff; if u_x < 0, forward diff.
        // Zero-gradient closure at domain boundary (no neighbour available).
        // Reference: Ferziger & Perić (2020) eq. 4.24 (donor-cell).
        //
        // dx = spacing between cell centres in x; derived from uniform grid layout.
        // For a uniform grid spanning [0, Lx] with Nx cells: centre(0,j).x = dx/2.
        const double dx = 2.0 * m_mesh->getCellCenter(0, 0).x();
        const double dy = 2.0 * m_mesh->getCellCenter(0, 0).y();

        for (int c = 0; c < N; ++c)
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

            bSrcX[c] = -(ux_c * dux_dx + uy_c * dux_dy) - (1.0 / rho) * gradP[c].x();
            bSrcY[c] = -(ux_c * duy_dx + uy_c * duy_dy) - (1.0 / rho) * gradP[c].y();
        }
    }

    // ── Build implicit viscous matrix ─────────────────────────────────────────
    // For component φ ∈ {ux, uy}:
    //   (V_P/dt + ν Σ_f a_f) φ*_P − ν a_f φ*_nb = (V_P/dt) φ^k_P + V_P bSrc_P
    //
    // a_f = A_face / |x_N − x_P|  for interior faces; boundary faces skipped (zero-gradient).
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
