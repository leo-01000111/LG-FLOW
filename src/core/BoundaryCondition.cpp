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

void BoundaryCondition::applyVelocity(Field<Eigen::Vector2d>& /*velocityField*/,
                                      const Mesh& /*mesh*/) const
{
    // STUB: Apply Dirichlet velocity at INLET/WALL faces and
    // zero-gradient at OUTLET/SYMMETRY. Will be implemented in Milestone 3.
}

void BoundaryCondition::applyPressure(Field<double>& /*pressureField*/,
                                      const Mesh& /*mesh*/) const
{
    // STUB: Apply fixed pressure at OUTLET, zero-gradient elsewhere.
    // Will be implemented in Milestone 3.
}
