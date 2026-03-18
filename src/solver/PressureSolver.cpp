#include "solver/PressureSolver.hpp"

PressureSolver::PressureSolver(const Mesh& mesh)
    : m_mesh(&mesh)
{
}

double PressureSolver::solve(Field<Eigen::Vector2d>& /*velocityField*/,
                             Field<double>&          /*pressureField*/,
                             double                  /*dt*/,
                             double                  /*rho*/,
                             double                  /*alphaP*/)
{
    // STUB: Pressure-correction solve (SIMPLE step 2) will be implemented in
    // Milestone 3. Returns zero residual for now.
    //
    // Implementation plan:
    //   1. Compute source term b_i = (rho/dt) * sum_f (u*_f · n_f) * A_f
    //   2. Assemble sparse coefficient matrix A from Laplacian stencil
    //   3. Solve A * p' = b using iterative solver
    //   4. Apply under-relaxation: p += alphaP * p'
    //   5. Correct velocity: u -= (dt/rho) * grad(p')
    return 0.0;
}
