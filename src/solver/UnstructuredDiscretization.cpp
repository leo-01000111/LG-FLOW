#include "solver/UnstructuredDiscretization.hpp"

#include <stdexcept>

// FVM Gauss-theorem operators on an unstructured 2D mesh.
// Interior faces: central differencing (2nd order).
// Boundary faces: zero-gradient closure (owner value extrapolated to face).
// Reference: Ferziger, Perić & Street — "Computational Methods for Fluid Dynamics",
//            4th ed. (2020), Chapters 4–6.

// ─────────────────────────────────────────────────────────────────────────────
// divergence  ∇·u ≈ (1/V) Σ_f (u_f · n_f) A_f   — Ferziger & Perić eq. 5.5
// ─────────────────────────────────────────────────────────────────────────────
UnstructuredField<double>
UnstructuredDiscretization::divergence(const UnstructuredField<Eigen::Vector2d>& field,
                                        const UnstructuredMesh&                   mesh)
{
    UnstructuredField<double> result(mesh, 0.0);
    const int nFaces = mesh.numFaces();

    for (int fi = 0; fi < nFaces; ++fi)
    {
        const MeshFace& f = mesh.face(FaceId{fi});
        const CellId P = f.owner;
        const CellId N = f.neighbour;

        const Eigen::Vector2d u_f = (toInt(N) < 0)
            ? field[P]
            : 0.5 * (field[P] + field[N]);

        const double flux = u_f.dot(f.normal) * f.area;

        result[P] += flux / mesh.cell(P).volume;

        if (toInt(N) >= 0)
            result[N] -= flux / mesh.cell(N).volume;
    }

    return result;
}

// ─────────────────────────────────────────────────────────────────────────────
// gradient  ∇φ ≈ (1/V) Σ_f φ_f n_f A_f   — Ferziger & Perić eq. 5.4
// ─────────────────────────────────────────────────────────────────────────────
UnstructuredField<Eigen::Vector2d>
UnstructuredDiscretization::gradient(const UnstructuredField<double>& field,
                                      const UnstructuredMesh&          mesh)
{
    UnstructuredField<Eigen::Vector2d> result(mesh, Eigen::Vector2d::Zero());
    const int nFaces = mesh.numFaces();

    for (int fi = 0; fi < nFaces; ++fi)
    {
        const MeshFace& f = mesh.face(FaceId{fi});
        const CellId P = f.owner;
        const CellId N = f.neighbour;

        const double phi_f = (toInt(N) < 0)
            ? field[P]
            : 0.5 * (field[P] + field[N]);

        const Eigen::Vector2d contrib = phi_f * f.normal * f.area;

        result[P] += contrib / mesh.cell(P).volume;

        if (toInt(N) >= 0)
            result[N] -= contrib / mesh.cell(N).volume;
    }

    return result;
}

// ─────────────────────────────────────────────────────────────────────────────
// divergenceRhieChow — Rhie-Chow (1983) face-interpolated divergence
// Prevents checkerboard pressure oscillations on collocated grids.
// Reference: Rhie & Chow (1983); Ferziger, Perić & Street (2020) Section 7.5.
// ─────────────────────────────────────────────────────────────────────────────
UnstructuredField<double>
UnstructuredDiscretization::divergenceRhieChow(
    const UnstructuredField<Eigen::Vector2d>& velocity,
    const UnstructuredField<double>&          pressure,
    const UnstructuredMesh&                   mesh,
    double                                    dt,
    double                                    rho)
{
    const UnstructuredField<Eigen::Vector2d> gradP = gradient(pressure, mesh);

    UnstructuredField<double> result(mesh, 0.0);
    const int nFaces = mesh.numFaces();
    const double dtOverRho = dt / rho;

    for (int fi = 0; fi < nFaces; ++fi)
    {
        const MeshFace& f = mesh.face(FaceId{fi});
        const CellId P = f.owner;
        const CellId N = f.neighbour;

        double faceFlux{};

        if (toInt(N) < 0)
        {
            // Boundary face: owner velocity (zero-gradient, no RC correction).
            faceFlux = velocity[P].dot(f.normal) * f.area;
        }
        else
        {
            // Interior face: apply Rhie-Chow correction.
            // Linearly interpolated normal velocity.
            const double uBar_dot_n = 0.5 * (velocity[P] + velocity[N]).dot(f.normal);

            // Compact face pressure gradient (two-point stencil) in normal direction.
            const Eigen::Vector2d& xP = mesh.cell(P).centre;
            const Eigen::Vector2d& xN = mesh.cell(N).centre;
            const double dist = (xN - xP).norm();

            const double compactGrad = (pressure[N] - pressure[P]) / dist;

            // Cell-centre gradient interpolated to face, projected on normal.
            const double interpGrad = 0.5 * (gradP[P] + gradP[N]).dot(f.normal);

            // RC correction removes checkerboard mode.
            const double rcFaceVel = uBar_dot_n - dtOverRho * (compactGrad - interpGrad);
            faceFlux = rcFaceVel * f.area;
        }

        result[P] += faceFlux / mesh.cell(P).volume;
        if (toInt(N) >= 0)
            result[N] -= faceFlux / mesh.cell(N).volume;
    }

    return result;
}
