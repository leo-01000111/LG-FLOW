#include "core/BoundaryCondition.hpp"

#include <algorithm>
#include <stdexcept>

namespace
{
/// Maps an outward boundary normal of the structured mesh to its patch name.
const char* patchNameForNormal(const Eigen::Vector2d& n)
{
    if (n.x() < -0.5) return "left";
    if (n.x() >  0.5) return "right";
    if (n.y() < -0.5) return "bottom";
    return "top";
}

/// Geometry shared by velocity and pressure boundary faces.
struct BoundaryFaceGeometry
{
    int             face;
    int             owner;
    Eigen::Vector2d normal;
    double          area;
    double          distance;
    Eigen::Vector2d centre;
};

/// Collects all boundary faces in ascending face order.
/// Owner centre → face centre distance is half the cell extent along the
/// normal, i.e. 0.5·V/A for a rectangular cell.
std::vector<BoundaryFaceGeometry> collectBoundaryFaces(const Mesh& mesh)
{
    std::vector<BoundaryFaceGeometry> faces;
    faces.reserve(static_cast<std::size_t>(2 * (mesh.Nx() + mesh.Ny())));

    const int Ny     = mesh.Ny();
    const int nFaces = mesh.numFaces();
    for (int f = 0; f < nFaces; ++f)
    {
        const auto [owner, neighbour] = mesh.getNeighbors(f);
        if (neighbour >= 0)
            continue;

        const Eigen::Vector2d n        = mesh.getFaceNormal(f);
        const double          area     = mesh.getFaceArea(f);
        const double          distance = 0.5 * mesh.getCellVolume(owner) / area;
        const Eigen::Vector2d centre   =
            mesh.getCellCenter(owner / Ny, owner % Ny) + distance * n;

        faces.push_back({f, owner, n, area, distance, centre});
    }
    return faces;
}
}  // namespace

void BoundaryCondition::addPatch(const std::string& name, const BoundaryPatch& patch)
{
    m_patches[name] = patch;
}

const BoundaryPatch& BoundaryCondition::getPatch(const std::string& name) const
{
    auto it = m_patches.find(name);
    if (it == m_patches.end())
        throw std::out_of_range("BoundaryCondition: unknown patch '" + name + "'");
    return it->second;
}

void BoundaryCondition::applyVelocity(Field<Eigen::Vector2d>& velocityField,
                                      const Mesh& mesh) const
{
    if (&velocityField.mesh() != &mesh)
        throw std::invalid_argument(
            "BoundaryCondition::applyVelocity: field mesh does not match supplied mesh");

    const int Nx = mesh.Nx();
    const int Ny = mesh.Ny();

    // Apply in order: left, right, bottom, top.
    // Top is applied last and wins at all top-row corners (lid-driven cavity convention).
    //
    // Velocity policy:
    //   WALL / INLET  → Dirichlet: set boundary cell to patch.value
    //   OUTLET / SYMMETRY → zero-gradient: leave boundary cell unchanged

    // Helper: apply Poiseuille profile u(y)=6·U_avg·y·(Ly-y)/Ly² on a column.
    // Ly is recovered from the first and last cell-centre y-coordinates.
    auto applyPoiseuille = [&](int col, double uAvg) {
        // Cell centres: y_j = (j+0.5)*dy, so Ly = y_0 + y_{Ny-1} = dy/2 + Ly-dy/2
        const double y0   = mesh.getCellCenter(col, 0).y();
        const double yLast= mesh.getCellCenter(col, Ny - 1).y();
        const double Ly   = y0 + yLast;          // = Ly_domain (exact for uniform grid)
        for (int j = 0; j < Ny; ++j) {
            const double y  = mesh.getCellCenter(col, j).y();
            const double ux = 6.0 * uAvg * y * (Ly - y) / (Ly * Ly);
            velocityField(col, j) = Eigen::Vector2d(ux, 0.0);
        }
    };

    if (auto it = m_patches.find("left"); it != m_patches.end()) {
        const BoundaryPatch& p = it->second;
        if (p.type == BoundaryType::WALL || p.type == BoundaryType::INLET) {
            for (int j = 0; j < Ny; ++j)
                velocityField(0, j) = p.value;
        } else if (p.type == BoundaryType::PARABOLIC_INLET) {
            applyPoiseuille(0, p.value.x());
        }
    }

    if (auto it = m_patches.find("right"); it != m_patches.end()) {
        const BoundaryPatch& p = it->second;
        if (p.type == BoundaryType::WALL || p.type == BoundaryType::INLET) {
            for (int j = 0; j < Ny; ++j)
                velocityField(Nx - 1, j) = p.value;
        } else if (p.type == BoundaryType::PARABOLIC_INLET) {
            applyPoiseuille(Nx - 1, p.value.x());
        }
    }

    if (auto it = m_patches.find("bottom"); it != m_patches.end()) {
        const BoundaryPatch& p = it->second;
        if (p.type == BoundaryType::WALL || p.type == BoundaryType::INLET) {
            for (int i = 0; i < Nx; ++i)
                velocityField(i, 0) = p.value;
        }
    }

    if (auto it = m_patches.find("top"); it != m_patches.end()) {
        const BoundaryPatch& p = it->second;
        if (p.type == BoundaryType::WALL || p.type == BoundaryType::INLET) {
            for (int i = 0; i < Nx; ++i)
                velocityField(i, Ny - 1) = p.value;
        }
    }
}

void BoundaryCondition::applyPressure(Field<double>& pressureField,
                                      const Mesh& mesh) const
{
    if (&pressureField.mesh() != &mesh)
        throw std::invalid_argument(
            "BoundaryCondition::applyPressure: field mesh does not match supplied mesh");

    const int Nx = mesh.Nx();
    const int Ny = mesh.Ny();

    // Apply in order: left, right, bottom, top.
    //
    // Pressure policy:
    //   OUTLET    → fixed reference: set boundary cell to 0.0
    //   INLET / WALL / SYMMETRY → zero-gradient: leave boundary cell unchanged

    if (auto it = m_patches.find("left"); it != m_patches.end()) {
        if (it->second.type == BoundaryType::OUTLET) {
            for (int j = 0; j < Ny; ++j)
                pressureField(0, j) = 0.0;
        }
    }

    if (auto it = m_patches.find("right"); it != m_patches.end()) {
        if (it->second.type == BoundaryType::OUTLET) {
            for (int j = 0; j < Ny; ++j)
                pressureField(Nx - 1, j) = 0.0;
        }
    }

    if (auto it = m_patches.find("bottom"); it != m_patches.end()) {
        if (it->second.type == BoundaryType::OUTLET) {
            for (int i = 0; i < Nx; ++i)
                pressureField(i, 0) = 0.0;
        }
    }

    if (auto it = m_patches.find("top"); it != m_patches.end()) {
        if (it->second.type == BoundaryType::OUTLET) {
            for (int i = 0; i < Nx; ++i)
                pressureField(i, Ny - 1) = 0.0;
        }
    }
}

std::vector<int> BoundaryCondition::collectCellsOfType(BoundaryType   type,
                                                        const Mesh&    mesh) const
{
    const int Nx = mesh.Nx();
    const int Ny = mesh.Ny();

    std::vector<int> cells;

    // Cell ID convention: c = i * Ny + j  (i = x-column, j = y-row)
    if (auto it = m_patches.find("left"); it != m_patches.end() && it->second.type == type)
        for (int j = 0; j < Ny; ++j)
            cells.push_back(0 * Ny + j);

    if (auto it = m_patches.find("right"); it != m_patches.end() && it->second.type == type)
        for (int j = 0; j < Ny; ++j)
            cells.push_back((Nx - 1) * Ny + j);

    if (auto it = m_patches.find("bottom"); it != m_patches.end() && it->second.type == type)
        for (int i = 0; i < Nx; ++i)
            cells.push_back(i * Ny + 0);

    if (auto it = m_patches.find("top"); it != m_patches.end() && it->second.type == type)
        for (int i = 0; i < Nx; ++i)
            cells.push_back(i * Ny + (Ny - 1));

    // Sort and deduplicate (corner cells may appear from two patches).
    std::sort(cells.begin(), cells.end());
    cells.erase(std::unique(cells.begin(), cells.end()), cells.end());
    return cells;
}

std::vector<VelocityBoundaryFace> BoundaryCondition::buildVelocityFaces(const Mesh& mesh) const
{
    const std::vector<BoundaryFaceGeometry> geometry = collectBoundaryFaces(mesh);

    // Channel height for PARABOLIC_INLET, recovered from the cell centres:
    // y_0 + y_{Ny-1} = dy/2 + (Ly - dy/2) = Ly on a uniform grid.
    const double Ly = mesh.getCellCenter(0, 0).y()
                    + mesh.getCellCenter(0, mesh.Ny() - 1).y();

    std::vector<VelocityBoundaryFace> faces;
    faces.reserve(geometry.size());

    for (const BoundaryFaceGeometry& g : geometry)
    {
        VelocityBoundaryFace vf{g.face, g.owner, g.normal, g.area, g.distance,
                                VelocityBCKind::ZERO_GRADIENT, Eigen::Vector2d::Zero()};

        const std::string name = patchNameForNormal(g.normal);
        if (auto it = m_patches.find(name); it != m_patches.end())
        {
            const BoundaryPatch& p = it->second;
            switch (p.type)
            {
            case BoundaryType::WALL:
            case BoundaryType::INLET:
                vf.kind  = VelocityBCKind::DIRICHLET;
                vf.value = p.value;
                break;
            case BoundaryType::PARABOLIC_INLET:
            {
                if (name == "top" || name == "bottom")
                    throw std::invalid_argument(
                        "BoundaryCondition::buildVelocityFaces: PARABOLIC_INLET is only "
                        "supported on 'left' and 'right' patches (got '" + name + "')");
                // Face-averaged Poiseuille profile, walls at y = 0 and y = Ly:
                //   (1/h)∫ 6U·y(Ly−y)/Ly² dy over [y_c − h/2, y_c + h/2]
                //     = 6U/Ly² · (y_c(Ly − y_c) − h²/12)
                // Averaging (not point-sampling) makes the inflow exactly U_avg·Ly.
                // Reference: White (2006) "Fluid Mechanics" §6-2.
                const double y = g.centre.y();
                const double h = g.area;
                vf.kind  = VelocityBCKind::DIRICHLET;
                vf.value = Eigen::Vector2d(
                    6.0 * p.value.x() * (y * (Ly - y) - h * h / 12.0) / (Ly * Ly), 0.0);
                break;
            }
            case BoundaryType::OUTLET:
            case BoundaryType::SYMMETRY:
                break;  // zero-gradient
            }
        }
        faces.push_back(vf);
    }
    return faces;
}

std::vector<PressureBoundaryFace> BoundaryCondition::buildPressureFaces(const Mesh& mesh) const
{
    const std::vector<BoundaryFaceGeometry> geometry = collectBoundaryFaces(mesh);

    std::vector<PressureBoundaryFace> faces;
    faces.reserve(geometry.size());

    for (const BoundaryFaceGeometry& g : geometry)
    {
        PressureBoundaryFace pf{g.face, g.owner, g.normal, g.area, g.distance,
                                PressureBCKind::ZERO_GRADIENT, 0.0};

        if (auto it = m_patches.find(patchNameForNormal(g.normal));
            it != m_patches.end() && it->second.type == BoundaryType::OUTLET)
        {
            pf.kind  = PressureBCKind::FIXED;
            pf.value = 0.0;
        }
        faces.push_back(pf);
    }
    return faces;
}
