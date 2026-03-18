#include "core/BoundaryCondition.hpp"

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

    if (auto it = m_patches.find("left"); it != m_patches.end()) {
        const BoundaryPatch& p = it->second;
        if (p.type == BoundaryType::WALL || p.type == BoundaryType::INLET) {
            for (int j = 0; j < Ny; ++j)
                velocityField(0, j) = p.value;
        }
    }

    if (auto it = m_patches.find("right"); it != m_patches.end()) {
        const BoundaryPatch& p = it->second;
        if (p.type == BoundaryType::WALL || p.type == BoundaryType::INLET) {
            for (int j = 0; j < Ny; ++j)
                velocityField(Nx - 1, j) = p.value;
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
