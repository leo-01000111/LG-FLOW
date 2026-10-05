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

/// Adds interior-face Gauss-gradient contributions (linear face interpolation).
void addInteriorGradient(const Field<double>& field, const Mesh& mesh,
                         Field<Eigen::Vector2d>& result)
{
    const int nFaces = mesh.numFaces();
    for (int f = 0; f < nFaces; ++f)
    {
        const auto [owner, neighbour] = mesh.getNeighbors(f);
        if (neighbour == -1)
            continue;

        const double          phi_f   = 0.5 * (field[owner] + field[neighbour]);
        const Eigen::Vector2d contrib = phi_f * mesh.getFaceNormal(f) * mesh.getFaceArea(f);

        result[owner]     += contrib / mesh.getCellVolume(owner);
        result[neighbour] -= contrib / mesh.getCellVolume(neighbour);
    }
}

/// Adds interior-face Rhie-Chow fluxes (Rhie & Chow 1983) to the divergence.
void addInteriorRhieChowFlux(const Field<Eigen::Vector2d>& velocity,
                             const Field<double>&          pressure,
                             const Field<Eigen::Vector2d>& gradP,
                             const Mesh&                   mesh,
                             double                        dtOverRho,
                             Field<double>&                result)
{
    const int nFaces = mesh.numFaces();
    for (int f = 0; f < nFaces; ++f)
    {
        const auto [owner, neighbour] = mesh.getNeighbors(f);
        if (neighbour == -1)
            continue;

        const Eigen::Vector2d n = mesh.getFaceNormal(f);  // owner → neighbour
        const double          A = mesh.getFaceArea(f);

        // Linearly interpolated normal face velocity:
        //   ū_f · n̂ = 0.5*(u_P + u_N) · n̂
        const double uBar_dot_n = 0.5 * (velocity[owner] + velocity[neighbour]).dot(n);

        // Compact face pressure gradient in the face-normal direction:
        //   (p_N − p_P) / |x_N − x_P|
        const double dist =
            (cellCentreById(neighbour, mesh) - cellCentreById(owner, mesh)).norm();
        const double compactGrad = (pressure[neighbour] - pressure[owner]) / dist;

        // Face-interpolated cell-centre gradient (projected on n̂):
        //   0.5*(∇p_P + ∇p_N) · n̂
        const double interpGrad = 0.5 * (gradP[owner] + gradP[neighbour]).dot(n);

        // Rhie-Chow correction: subtract (dt/ρ) * (compact − interp).
        // The compact gradient sees the checkerboard; the interpolated gradient
        // does not.  The difference damps decoupled pressure modes.
        const double faceFlux = (uBar_dot_n - dtOverRho * (compactGrad - interpGrad)) * A;

        result[owner]     += faceFlux / mesh.getCellVolume(owner);
        result[neighbour] -= faceFlux / mesh.getCellVolume(neighbour);
    }
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
// gradient with boundary-face closure — Ferziger & Perić eq. 5.4
// Central interpolation on interior faces; BC-prescribed values at boundary faces.
// ─────────────────────────────────────────────────────────────────────────────
Field<Eigen::Vector2d>
Discretization::gradient(const Field<double>&                   pressure,
                         const Mesh&                            mesh,
                         std::span<const PressureBoundaryFace>  boundaryFaces)
{
    if (&pressure.mesh() != &mesh)
        throw std::invalid_argument(
            "Discretization::gradient: field is not on the supplied mesh");

    Field<Eigen::Vector2d> result(mesh, Eigen::Vector2d::Zero());
    addInteriorGradient(pressure, mesh, result);

    for (const PressureBoundaryFace& bf : boundaryFaces)
    {
        const double phi_f = (bf.kind == PressureBCKind::FIXED) ? bf.value
                                                                : pressure[bf.owner];
        result[bf.owner] += phi_f * bf.normal * bf.area / mesh.getCellVolume(bf.owner);
    }
    return result;
}

Field<Eigen::Vector2d>
Discretization::gradient(const Field<double>&                   component,
                         const Mesh&                            mesh,
                         std::span<const VelocityBoundaryFace>  boundaryFaces,
                         int                                    axis)
{
    if (&component.mesh() != &mesh)
        throw std::invalid_argument(
            "Discretization::gradient: field is not on the supplied mesh");
    if (axis != 0 && axis != 1)
        throw std::invalid_argument(
            "Discretization::gradient: axis must be 0 or 1 (got " + std::to_string(axis) + ")");

    Field<Eigen::Vector2d> result(mesh, Eigen::Vector2d::Zero());
    addInteriorGradient(component, mesh, result);

    for (const VelocityBoundaryFace& bf : boundaryFaces)
    {
        const double phi_f = (bf.kind == VelocityBCKind::DIRICHLET) ? bf.value[axis]
                                                                    : component[bf.owner];
        result[bf.owner] += phi_f * bf.normal * bf.area / mesh.getCellVolume(bf.owner);
    }
    return result;
}

// ─────────────────────────────────────────────────────────────────────────────
// divergenceRhieChow — Rhie-Chow (1983) face-interpolated divergence
// Prevents checkerboard pressure oscillations on collocated grids by replacing
// the linearly-interpolated face velocity with a compact-gradient-corrected value.
// Reference: Rhie & Chow (1983); Ferziger, Perić & Street (2020) Section 7.5.
// ─────────────────────────────────────────────────────────────────────────────
Field<double>
Discretization::divergenceRhieChow(const Field<Eigen::Vector2d>& velocity,
                                    const Field<double>&          pressure,
                                    const Mesh&                   mesh,
                                    double                        dt,
                                    double                        rho)
{
    if (&velocity.mesh() != &mesh)
        throw std::invalid_argument(
            "Discretization::divergenceRhieChow: velocity field is not on the supplied mesh");
    if (&pressure.mesh() != &mesh)
        throw std::invalid_argument(
            "Discretization::divergenceRhieChow: pressure field is not on the supplied mesh");

    // Cell-centre pressure gradients for face interpolation.
    // Used in the RC correction: 0.5*(∇p_P + ∇p_N) · n̂_f.
    const Field<Eigen::Vector2d> gradP = gradient(pressure, mesh);

    Field<double> result(mesh, 0.0);
    addInteriorRhieChowFlux(velocity, pressure, gradP, mesh, dt / rho, result);

    // Boundary faces: use owner velocity (zero-gradient; same as divergence()).
    const int nFaces = mesh.numFaces();
    for (int f = 0; f < nFaces; ++f)
    {
        const auto [owner, neighbour] = mesh.getNeighbors(f);
        if (neighbour != -1)
            continue;
        const double faceFlux =
            velocity[owner].dot(mesh.getFaceNormal(f)) * mesh.getFaceArea(f);
        result[owner] += faceFlux / mesh.getCellVolume(owner);
    }

    return result;
}

Field<double>
Discretization::divergenceRhieChow(const Field<Eigen::Vector2d>&          velocity,
                                    const Field<double>&                   pressure,
                                    const Mesh&                            mesh,
                                    double                                 dt,
                                    double                                 rho,
                                    std::span<const VelocityBoundaryFace>  velocityFaces,
                                    std::span<const PressureBoundaryFace>  pressureFaces)
{
    if (&velocity.mesh() != &mesh)
        throw std::invalid_argument(
            "Discretization::divergenceRhieChow: velocity field is not on the supplied mesh");
    if (&pressure.mesh() != &mesh)
        throw std::invalid_argument(
            "Discretization::divergenceRhieChow: pressure field is not on the supplied mesh");
    if (velocityFaces.size() != pressureFaces.size())
        throw std::invalid_argument(
            "Discretization::divergenceRhieChow: velocity and pressure face lists differ in size");

    const Field<Eigen::Vector2d> gradP = gradient(pressure, mesh, pressureFaces);
    const double dtOverRho = dt / rho;

    Field<double> result(mesh, 0.0);
    addInteriorRhieChowFlux(velocity, pressure, gradP, mesh, dtOverRho, result);

    for (std::size_t k = 0; k < velocityFaces.size(); ++k)
    {
        const VelocityBoundaryFace& vf = velocityFaces[k];
        const PressureBoundaryFace& pf = pressureFaces[k];
        if (vf.face != pf.face)
            throw std::invalid_argument(
                "Discretization::divergenceRhieChow: face lists are not in the same order");

        const int owner = vf.owner;
        double faceVel;  // normal face velocity û_f · n̂

        if (vf.kind == VelocityBCKind::DIRICHLET)
        {
            // Prescribed velocity: the boundary flux is known exactly.
            faceVel = vf.value.dot(vf.normal);
        }
        else if (pf.kind == PressureBCKind::FIXED)
        {
            // Outlet: Rhie-Chow between owner and face with the fixed face
            // pressure, consistent with the p' = 0 face condition in PressureSolver.
            const double compactGrad = (pf.value - pressure[owner]) / pf.distance;
            const double cellGrad    = gradP[owner].dot(vf.normal);
            faceVel = velocity[owner].dot(vf.normal) - dtOverRho * (compactGrad - cellGrad);
        }
        else
        {
            faceVel = velocity[owner].dot(vf.normal);
        }

        result[owner] += faceVel * vf.area / mesh.getCellVolume(owner);
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
