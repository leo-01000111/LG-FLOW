#include "solver/PressureSolver.hpp"
#include "solver/Discretization.hpp"

#include <Eigen/Sparse>
#include <Eigen/SparseLU>

#include <cmath>
#include <stdexcept>
#include <vector>

PressureSolver::PressureSolver(const Mesh& mesh)
    : m_mesh(&mesh)
{
}

double PressureSolver::solve(Field<Eigen::Vector2d>& velocityField,
                             Field<double>&          pressureField,
                             double                  dt,
                             double                  rho,
                             double                  alphaP)
{
    // --- Input validation ---
    if (dt <= 0.0)
        throw std::invalid_argument("PressureSolver::solve: dt must be > 0");
    if (rho <= 0.0)
        throw std::invalid_argument("PressureSolver::solve: rho must be > 0");
    if (alphaP <= 0.0 || alphaP > 1.0)
        throw std::invalid_argument("PressureSolver::solve: alphaP must be in (0, 1]");
    if (&velocityField.mesh() != m_mesh)
        throw std::invalid_argument("PressureSolver::solve: velocityField mesh mismatch");
    if (&pressureField.mesh() != m_mesh)
        throw std::invalid_argument("PressureSolver::solve: pressureField mesh mismatch");

    const int N  = m_mesh->numCells();
    const int Nx = m_mesh->Nx();
    const int Ny = m_mesh->Ny();

    // --- Step 1: Build RHS source term ---
    // b_i = (rho/dt) * divergence(u*)_i * V_i
    // divergence() returns (1/V) Σ_f (u_f · n_f) A_f; multiply back by V.
    // Reference: Patankar (1980) eq. 6.28
    Field<double> divU = Discretization::divergence(velocityField, *m_mesh);

    Eigen::VectorXd b(N);
    for (int cellId = 0; cellId < N; ++cellId)
    {
        const double V = m_mesh->getCellVolume(cellId);
        // A is positive-definite (a_f > 0 diagonal dominant), so
        // (A p')[o] = -Laplacian(p')[o] * V_o.  To get ∇²p' = (ρ/dt)·div(u*),
        // the RHS must carry a negative sign.
        b[cellId] = -(rho / dt) * divU[cellId] * V;
    }

    // --- Step 2: Assemble sparse Poisson-like pressure-correction matrix ---
    // 5-point stencil, cell-centred, face-normal diffusion form.
    // For each interior face f with owner o and neighbour n:
    //   a_f = faceArea / distance(centre_o, centre_n)
    //   A[o,o] += a_f;  A[n,n] += a_f
    //   A[o,n] -= a_f;  A[n,o] -= a_f
    // Boundary faces (n == -1): zero-normal-gradient → no off-diagonal coupling.
    // Reference: Ferziger, Peric & Street, 4th ed., eq. 7.23
    using Triplet = Eigen::Triplet<double>;
    std::vector<Triplet> triplets;
    triplets.reserve(static_cast<std::size_t>(5 * N));

    std::vector<double> diag(static_cast<std::size_t>(N), 0.0);

    const int nFaces = m_mesh->numFaces();
    for (int f = 0; f < nFaces; ++f)
    {
        const auto [o, n] = m_mesh->getNeighbors(f);
        if (n < 0)
            continue;  // boundary face: zero-gradient, no coupling

        const int io  = o / Ny;
        const int jo  = o % Ny;
        const int in_ = n / Ny;
        const int jn  = n % Ny;

        const Eigen::Vector2d xo = m_mesh->getCellCenter(io, jo);
        const Eigen::Vector2d xn = m_mesh->getCellCenter(in_, jn);

        const double area = m_mesh->getFaceArea(f);
        const double dist = (xn - xo).norm();
        const double aF   = area / dist;

        diag[static_cast<std::size_t>(o)] += aF;
        diag[static_cast<std::size_t>(n)] += aF;

        // Skip row-0 off-diagonal: reference cell row is overridden below.
        // Column-0 coupling in other rows is harmless (x[0] is forced to 0).
        if (o != 0)
            triplets.emplace_back(o, n, -aF);
        if (n != 0)
            triplets.emplace_back(n, o, -aF);
    }

    // Diagonal entries — skip cell 0 (overridden by identity row below)
    for (int cellId = 1; cellId < N; ++cellId)
        triplets.emplace_back(cellId, cellId, diag[static_cast<std::size_t>(cellId)]);

    // --- Step 3: Fix pure-Neumann singularity ---
    // Row 0 → identity; b[0] = 0 forces p'[0] = 0 as the reference cell.
    // BiCGSTAB tolerates the resulting asymmetric row.
    triplets.emplace_back(0, 0, 1.0);
    b[0] = 0.0;

    Eigen::SparseMatrix<double> A(N, N);
    A.setFromTriplets(triplets.begin(), triplets.end());
    A.makeCompressed();

    // --- Step 4: Solve A * x = b (SparseLU direct solver) ---
    // SparseLU is used instead of BiCGSTAB because the row-0 identity fix
    // creates an asymmetric matrix that is difficult for iterative solvers
    // at large mesh sizes (64x64 and above). SparseLU gives exact results
    // for the sparse 5-point Laplacian stencil without iteration count risk.
    // Reference: Eigen documentation — SparseLU.
    Eigen::SparseLU<Eigen::SparseMatrix<double>> directSolver;
    directSolver.analyzePattern(A);
    directSolver.factorize(A);

    if (directSolver.info() != Eigen::Success)
        throw std::runtime_error(
            "PressureSolver::solve: SparseLU factorization failed (singular matrix?)");

    const Eigen::VectorXd x = directSolver.solve(b);

    if (directSolver.info() != Eigen::Success)
        throw std::runtime_error(
            "PressureSolver::solve: SparseLU solve step failed");

    // --- Step 5: Verify finite solution ---
    for (int i = 0; i < N; ++i)
    {
        if (!std::isfinite(x[i]))
            throw std::runtime_error(
                "PressureSolver::solve: non-finite pressure correction at cell " +
                std::to_string(i));
    }

    // --- Step 6: Build p' field from solution vector ---
    Field<double> pPrime(*m_mesh, 0.0);
    for (int cellId = 0; cellId < N; ++cellId)
        pPrime[cellId] = x[cellId];

    // --- Step 7: Correct pressure ---
    // p ← p + αp * p'   (Patankar 1980, eq. 6.30)
    for (int cellId = 0; cellId < N; ++cellId)
        pressureField[cellId] += alphaP * pPrime[cellId];

    // --- Step 8: Correct velocity (stencil-consistent structured-grid gradient) ---
    // u ← u* − (dt/ρ) ∇p'
    //
    // Replaces the generic Discretization::gradient (which iterates over the
    // face list) with an equivalent inline computation on the structured grid.
    // This avoids repeated face-list traversal and makes the stencil explicit:
    //
    //   Interior cells (both neighbours exist): central difference
    //     ∂p'/∂x ≈ (p'[i+1,j] − p'[i-1,j]) / (2·dx)
    //   Boundary cells (one neighbour missing): zero-gradient closure
    //     face value = cell value → same half-difference as Gauss theorem gives
    //     ∂p'/∂x|_{i=0}   ≈ (p'[1,j]   − p'[0,j]) / (2·dx)
    //     ∂p'/∂x|_{i=Nx-1}≈ (p'[Nx-1,j]− p'[Nx-2,j]) / (2·dx)
    //
    // This is algebraically identical to Discretization::gradient but avoids the
    // general face-iteration overhead and keeps the stencil readable.
    // Reference: Patankar (1980) eq. 6.31; Ferziger & Peric (2020) eq. 4.22.
    {
        const double coeff = dt / rho;
        // Derive uniform cell spacing from the first cell centre.
        // For a uniform grid spanning [0, Lx]: centre(0,0) = (dx/2, dy/2).
        const double dx = 2.0 * m_mesh->getCellCenter(0, 0).x();
        const double dy = 2.0 * m_mesh->getCellCenter(0, 0).y();

        for (int cellId = 0; cellId < N; ++cellId)
        {
            const int i = cellId / Ny;
            const int j = cellId % Ny;

            // x-gradient: central for interior cells, zero-gradient closure
            // at boundaries (boundary face value = cell centre value).
            double gradX = 0.0;
            if (Nx > 1)
            {
                // Boundary: face value p'_boundary = p'_cell → half-diff
                const double pL = (i > 0)      ? pPrime[(i - 1) * Ny + j] : pPrime[cellId];
                const double pR = (i < Nx - 1) ? pPrime[(i + 1) * Ny + j] : pPrime[cellId];
                gradX = (pR - pL) / (2.0 * dx);
            }

            // y-gradient: same pattern, y-direction.
            double gradY = 0.0;
            if (Ny > 1)
            {
                const double pB = (j > 0)      ? pPrime[i * Ny + j - 1] : pPrime[cellId];
                const double pT = (j < Ny - 1) ? pPrime[i * Ny + j + 1] : pPrime[cellId];
                gradY = (pT - pB) / (2.0 * dy);
            }

            velocityField[cellId].x() -= coeff * gradX;
            velocityField[cellId].y() -= coeff * gradY;
        }
    }

    // --- Step 9: Residual ||A x − b||_2 ---
    const Eigen::VectorXd residualVec = A * x - b;
    const double residual = residualVec.norm();

    if (!std::isfinite(residual))
        throw std::runtime_error(
            "PressureSolver::solve: non-finite residual norm");

    return residual;
}
