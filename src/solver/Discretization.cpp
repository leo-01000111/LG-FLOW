#include "solver/Discretization.hpp"

#include <stdexcept>
#include <string>

// FVM Gauss-theorem operators on a structured uniform Cartesian mesh.
// All face loops iterate over the global face list provided by Mesh.
// Interior faces use central differencing (2nd order).
// Boundary faces use zero-gradient closure (owner value extrapolated to face).
// Reference: Ferziger, Perić & Street — "Computational Methods for Fluid Dynamics",
//            4th ed. (2020), Chapters 4–6.

namespace
{
/// Returns the cell-centre coordinates for a linear cell index.
inline Eigen::Vector2d cellCentreById(int cellId, const Mesh& mesh)
{
    const int i = cellId / mesh.Ny();
    const int j = cellId % mesh.Ny();
    return mesh.getCellCenter(i, j);
}
}  // namespace

// ─────────────────────────────────────────────────────────────────────────────
// divergence  ∇·u ≈ (1/V) Σ_f (u_f · n_f) A_f   — Ferziger & Perić eq. 5.5
// Central differencing on interior faces; zero-gradient at boundary faces.
// ─────────────────────────────────────────────────────────────────────────────
Field<double>
Discretization::divergence(const Field<Eigen::Vector2d>& field, const Mesh& mesh)
{
    if (&field.mesh() != &mesh)
        throw std::invalid_argument(
            "Discretization::divergence: field is not on the supplied mesh");

    Field<double> result(mesh, 0.0);
    const int nFaces = mesh.numFaces();

    for (int f = 0; f < nFaces; ++f)
    {
        const auto [owner, neighbour] = mesh.getNeighbors(f);
        const Eigen::Vector2d n = mesh.getFaceNormal(f);
        const double A = mesh.getFaceArea(f);

        // Face velocity: zero-gradient at boundary, central for interior
        const Eigen::Vector2d u_f = (neighbour == -1)
            ? field[owner]
            : 0.5 * (field[owner] + field[neighbour]);

        const double flux = u_f.dot(n) * A;

        result[owner] += flux / mesh.getCellVolume(owner);

        if (neighbour != -1)
            result[neighbour] -= flux / mesh.getCellVolume(neighbour);
    }

    return result;
}

// ─────────────────────────────────────────────────────────────────────────────
// gradient  ∇φ ≈ (1/V) Σ_f φ_f n_f A_f   — Ferziger & Perić eq. 5.4
// Central differencing on interior faces; zero-gradient at boundary faces.
// ─────────────────────────────────────────────────────────────────────────────
Field<Eigen::Vector2d>
Discretization::gradient(const Field<double>& field, const Mesh& mesh)
{
    if (&field.mesh() != &mesh)
        throw std::invalid_argument(
            "Discretization::gradient: field is not on the supplied mesh");

    Field<Eigen::Vector2d> result(mesh, Eigen::Vector2d::Zero());
    const int nFaces = mesh.numFaces();

    for (int f = 0; f < nFaces; ++f)
    {
        const auto [owner, neighbour] = mesh.getNeighbors(f);
        const Eigen::Vector2d n = mesh.getFaceNormal(f);
        const double A = mesh.getFaceArea(f);

        // Face scalar value: zero-gradient at boundary, central for interior
        const double phi_f = (neighbour == -1)
            ? field[owner]
            : 0.5 * (field[owner] + field[neighbour]);

        // Gauss theorem: outward-normal contribution to owner is +n, to neighbour is -n
        const Eigen::Vector2d contrib = phi_f * n * A;

        result[owner] += contrib / mesh.getCellVolume(owner);

        if (neighbour != -1)
            result[neighbour] -= contrib / mesh.getCellVolume(neighbour);
    }

    return result;
}

// ─────────────────────────────────────────────────────────────────────────────
// laplacian  ∇²φ ≈ (1/V) Σ_f (∂φ/∂n_f) A_f   — Ferziger & Perić eq. 4.22
// Compact two-point stencil on interior faces; zero contribution at boundaries.
// ─────────────────────────────────────────────────────────────────────────────
Field<double>
Discretization::laplacian(const Field<double>& field, const Mesh& mesh)
{
    if (&field.mesh() != &mesh)
        throw std::invalid_argument(
            "Discretization::laplacian: field is not on the supplied mesh");

    Field<double> result(mesh, 0.0);
    const int nFaces = mesh.numFaces();

    for (int f = 0; f < nFaces; ++f)
    {
        const auto [owner, neighbour] = mesh.getNeighbors(f);
        if (neighbour == -1) continue;  // zero normal gradient at boundary faces

        const Eigen::Vector2d n = mesh.getFaceNormal(f);
        const double A = mesh.getFaceArea(f);

        // Normal distance from owner centre to neighbour centre along face normal
        const Eigen::Vector2d x_o = cellCentreById(owner, mesh);
        const Eigen::Vector2d x_n = cellCentreById(neighbour, mesh);
        const double d = (x_n - x_o).dot(n);

        if (d <= 0.0)
            throw std::runtime_error(
                "Discretization::laplacian: non-positive face distance d=" +
                std::to_string(d) + " at face " + std::to_string(f));

        // Normal gradient ≈ (φ_n − φ_o) / d  (central differencing)
        const double grad_n = (field[neighbour] - field[owner]) / d;
        const double flux = grad_n * A;

        result[owner] += flux / mesh.getCellVolume(owner);
        result[neighbour] -= flux / mesh.getCellVolume(neighbour);
    }

    return result;
}
