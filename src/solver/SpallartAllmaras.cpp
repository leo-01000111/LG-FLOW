#include "solver/SpallartAllmaras.hpp"

#include "solver/Discretization.hpp"

#include <Eigen/Sparse>
#include <Eigen/SparseLU>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

// ── SA model constants ────────────────────────────────────────────────────────
// Reference: Spalart & Allmaras (1994); Allmaras, Johnson & Spalart (2012).

namespace
{
constexpr double SA_CB1   = 0.1355;
constexpr double SA_CB2   = 0.622;
constexpr double SA_SIGMA = 2.0 / 3.0;
constexpr double SA_KAPPA = 0.41;
constexpr double SA_CW1   = SA_CB1 / (SA_KAPPA * SA_KAPPA) + (1.0 + SA_CB2) / SA_SIGMA;
constexpr double SA_CW2   = 0.3;
constexpr double SA_CW3   = 2.0;
constexpr double SA_CV1   = 7.1;
constexpr double SA_R_LIM = 10.0;  // clamp for r function (Allmaras 2012 §2.3)
} // namespace

// ── Constructor ───────────────────────────────────────────────────────────────

SpallartAllmaras::SpallartAllmaras(const Mesh& mesh)
    : m_mesh(&mesh)
{
}

// ── Closure functions ─────────────────────────────────────────────────────────

double SpallartAllmaras::chi(double nuTilde, double nu)
{
    return nuTilde / nu;
}

double SpallartAllmaras::fv1(double chiVal)
{
    const double chi3 = chiVal * chiVal * chiVal;
    const double cv13 = SA_CV1 * SA_CV1 * SA_CV1;
    return chi3 / (chi3 + cv13);
}

double SpallartAllmaras::fv2(double chiVal, double fv1Val)
{
    return 1.0 - chiVal / (1.0 + chiVal * fv1Val);
}

double SpallartAllmaras::sTilde(double omega, double d, double nu, double nuTilde)
{
    // Avoid division by zero near wall; return plain vorticity when d is negligible.
    if (d < 1.0e-14)
        return omega;

    const double nuT_pos = std::max(nuTilde, 0.0);
    const double chiVal  = chi(nuT_pos, nu);
    const double fv1Val  = fv1(chiVal);
    const double fv2Val  = fv2(chiVal, fv1Val);
    return omega + (nuT_pos / (SA_KAPPA * SA_KAPPA * d * d)) * fv2Val;
}

double SpallartAllmaras::rFunc(double nuTilde, double sTildeVal, double d)
{
    if (sTildeVal < 1.0e-14 || d < 1.0e-14)
        return SA_R_LIM;

    const double r = nuTilde / (sTildeVal * SA_KAPPA * SA_KAPPA * d * d);
    return std::min(r, SA_R_LIM);
}

double SpallartAllmaras::gFunc(double rVal)
{
    return rVal + SA_CW2 * (std::pow(rVal, 6.0) - rVal);
}

double SpallartAllmaras::fw(double gVal)
{
    const double cw3_6 = std::pow(SA_CW3, 6.0);
    const double g6    = std::pow(gVal,   6.0);
    return gVal * std::pow((1.0 + cw3_6) / (g6 + cw3_6), 1.0 / 6.0);
}

// ── Transport equation solve ──────────────────────────────────────────────────

void SpallartAllmaras::solve(Field<double>&                nuTilde,
                              const Field<Eigen::Vector2d>& velocity,
                              const Field<double>&          wallDist,
                              double                        nu,
                              double                        dt,
                              double                        alphaNu)
{
    const int N  = m_mesh->numCells();
    const int Ny = m_mesh->Ny();

    // Save previous iterate for under-relaxation reference.
    // Pre-allocate: copy before building the system, not inside any loop.
    std::vector<double> nuTildeOld(static_cast<std::size_t>(N));
    for (int c = 0; c < N; ++c)
        nuTildeOld[static_cast<std::size_t>(c)] = nuTilde[c];

    // ── Step 1: Vorticity magnitude Ω = |∂vy/∂x - ∂vx/∂y| ───────────────────
    // Decompose velocity into scalar components for Discretization::gradient.
    Field<double> ux(*m_mesh, 0.0);
    Field<double> uy(*m_mesh, 0.0);
    for (int c = 0; c < N; ++c)
    {
        ux[c] = velocity[c].x();
        uy[c] = velocity[c].y();
    }

    // Central-differencing gradient (Gauss theorem, Ferziger & Perić eq. 5.4).
    const Field<Eigen::Vector2d> gradUx = Discretization::gradient(ux, *m_mesh);
    const Field<Eigen::Vector2d> gradUy = Discretization::gradient(uy, *m_mesh);

    // ── Step 2: ∇ν̃ gradient ──────────────────────────────────────────────────
    const Field<Eigen::Vector2d> gradNuT = Discretization::gradient(nuTilde, *m_mesh);

    // ── Steps 3 & 4: Assemble matrix ─────────────────────────────────────────
    // Same sparse-triplet pattern as MomentumSolver (Ferziger & Perić §7.4).
    using Triplet = Eigen::Triplet<double>;
    std::vector<Triplet> triplets;
    triplets.reserve(static_cast<std::size_t>(5 * N));

    std::vector<double> diag(static_cast<std::size_t>(N), 0.0);
    Eigen::VectorXd     rhs(N);

    for (int c = 0; c < N; ++c)
    {
        const double V     = m_mesh->getCellVolume(c);
        const double d     = wallDist[c];
        const double nuT_c = std::max(nuTilde[c], 0.0);

        // Closure values at cell c
        const double omegaC    = std::abs(gradUy[c].x() - gradUx[c].y());
        const double sTildeVal = sTilde(omegaC, d, nu, nuT_c);
        const double rVal      = rFunc(nuT_c, sTildeVal, d);
        const double gVal      = gFunc(rVal);
        const double fwVal     = fw(gVal);

        // Transient (Euler implicit): V/dt on diagonal, (V/dt)*ν̃^k on RHS
        diag[static_cast<std::size_t>(c)] = V / dt;
        rhs[c]                            = (V / dt) * nuTildeOld[static_cast<std::size_t>(c)];

        // Explicit production: cb1 · S̃ · ν̃  (SA eq. 1, Spalart & Allmaras 1994)
        rhs[c] += V * SA_CB1 * sTildeVal * nuT_c;

        // Explicit cross-diffusion: (cb2/σ) · |∇ν̃|²  (SA eq. 1)
        rhs[c] += V * (SA_CB2 / SA_SIGMA) * gradNuT[c].squaredNorm();

        // Linearized destruction: cw1 · fw · ν̃^k / d²  → implicit (added to diagonal)
        // Linearization about ν̃^k prevents negative solutions amplifying destruction.
        const double dSq = std::max(d * d, 1.0e-14);
        diag[static_cast<std::size_t>(c)] += V * SA_CW1 * fwVal * nuT_c / dSq;
    }

    // ── Step 5: Viscous face coupling (implicit central) ─────────────────────
    // Diffusion coefficient at face f: (ν + ν̃_f) / σ
    // Face average: ν̃_f = 0.5 * (ν̃_owner + ν̃_neighbour)
    const int nFaces = m_mesh->numFaces();
    for (int f = 0; f < nFaces; ++f)
    {
        const auto [o, n] = m_mesh->getNeighbors(f);
        if (n < 0)
            continue;  // boundary face: zero-gradient, no coupling

        // Arithmetic mean of effective viscosity at the face: (ν + ν̃)/σ
        const double nuEff_f = 0.5 * ((nu + nuTilde[o]) + (nu + nuTilde[n]));

        const int io  = o / Ny;
        const int jo  = o % Ny;
        const int in_ = n / Ny;
        const int jn  = n % Ny;

        const Eigen::Vector2d xo   = m_mesh->getCellCenter(io,  jo);
        const Eigen::Vector2d xn   = m_mesh->getCellCenter(in_, jn);
        const double          area = m_mesh->getFaceArea(f);
        const double          dist = (xn - xo).norm();
        const double          aF   = (nuEff_f / SA_SIGMA) * area / dist;

        diag[static_cast<std::size_t>(o)] += aF;
        diag[static_cast<std::size_t>(n)] += aF;

        triplets.emplace_back(o, n, -aF);
        triplets.emplace_back(n, o, -aF);
    }

    // ── Step 6: Explicit upwind convection for ν̃ ─────────────────────────────
    // First-order upwind: face ν̃ = owner ν̃ if flux > 0, neighbour ν̃ if flux < 0.
    // The outward normal from owner defines the flux sign convention.
    for (int f = 0; f < nFaces; ++f)
    {
        const auto [o, n] = m_mesh->getNeighbors(f);
        if (n < 0)
            continue;

        const Eigen::Vector2d uFace    = 0.5 * (velocity[o] + velocity[n]);
        const Eigen::Vector2d faceNorm = m_mesh->getFaceNormal(f);
        const double          area     = m_mesh->getFaceArea(f);
        // Volumetric flux (positive = out of owner cell)
        const double flux = uFace.dot(faceNorm) * area;

        const double convNuT = (flux > 0.0) ? nuTilde[o] : nuTilde[n];

        // RHS[owner] -= conv_nuT * flux (outgoing flux removes ν̃ from owner)
        // RHS[neighbour] += conv_nuT * flux (incoming flux adds ν̃ to neighbour)
        rhs[o] -= convNuT * flux;
        rhs[n] += convNuT * flux;
    }

    // ── Step 7: Under-relaxation (Patankar 1980, eq. 6.36) ───────────────────
    // Replace a_P with a_P/α and add (1-α)/α · a_P · ν̃^k to RHS.
    if (alphaNu < 1.0)
    {
        for (int c = 0; c < N; ++c)
        {
            const double aP    = diag[static_cast<std::size_t>(c)];
            const double relax = (1.0 - alphaNu) / alphaNu;
            rhs[c] += relax * aP * nuTildeOld[static_cast<std::size_t>(c)];
            diag[static_cast<std::size_t>(c)] = aP / alphaNu;
        }
    }

    // ── Step 8: Fill diagonal and solve ──────────────────────────────────────
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
            "SpallartAllmaras::solve: SparseLU factorization failed");

    const Eigen::VectorXd sol = linSolver.solve(rhs);
    if (linSolver.info() != Eigen::Success)
        throw std::runtime_error(
            "SpallartAllmaras::solve: SparseLU solve failed");

    // ── Step 9: Write solution, clamp negatives, check finite ─────────────────
    for (int c = 0; c < N; ++c)
    {
        if (!std::isfinite(sol[c]))
            throw std::runtime_error(
                "SpallartAllmaras::solve: non-finite nu_tilde at cell "
                + std::to_string(c));

        // Allmaras 2012 §2.4: negative ν̃ is non-physical; clamp to zero.
        nuTilde[c] = std::max(sol[c], 0.0);
    }
}

// ── computeNuT ────────────────────────────────────────────────────────────────

Field<double> SpallartAllmaras::computeNuT(const Field<double>& nuTilde,
                                            double               nu,
                                            const Mesh&          mesh)
{
    const int     N = mesh.numCells();
    Field<double> nuT(mesh, 0.0);

    for (int c = 0; c < N; ++c)
    {
        if (nuTilde[c] > 0.0)
        {
            const double chiVal = chi(nuTilde[c], nu);
            nuT[c] = nuTilde[c] * fv1(chiVal);
        }
        // nuTilde <= 0: nuT stays 0 (Allmaras 2012 negative-branch fix)
    }

    return nuT;
}
