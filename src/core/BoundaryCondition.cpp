#include "core/BoundaryCondition.hpp"

#include <algorithm>
#include <stdexcept>

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
