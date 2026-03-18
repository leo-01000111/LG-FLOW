#include "solver/Discretization.hpp"

// STUB: All operators return zero-initialised fields.
// Full FVM discretization (central differencing, Gauss theorem) will be
// implemented in Milestone 3.

Field<double>
Discretization::divergence(const Field<Eigen::Vector2d>& field, const Mesh& mesh)
{
    // STUB: ∇·u via Gauss theorem — Ferziger & Perić eq. 5.5
    return Field<double>(mesh, 0.0);
}

Field<Eigen::Vector2d>
Discretization::gradient(const Field<double>& field, const Mesh& mesh)
{
    // STUB: ∇φ via Gauss theorem — Ferziger & Perić eq. 5.4
    return Field<Eigen::Vector2d>(mesh, Eigen::Vector2d::Zero());
}

Field<double>
Discretization::laplacian(const Field<double>& field, const Mesh& mesh)
{
    // STUB: ∇²φ via compact two-point stencil — Ferziger & Perić eq. 4.22
    return Field<double>(mesh, 0.0);
}
