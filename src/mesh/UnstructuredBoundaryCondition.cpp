#include "mesh/UnstructuredBoundaryCondition.hpp"

#include <algorithm>
#include <stdexcept>

void UnstructuredBoundaryCondition::addPatch(const std::string& name,
                                              const BoundaryPatch& patch)
{
    m_patches[name] = patch;
}

const BoundaryPatch& UnstructuredBoundaryCondition::getPatch(const std::string& name) const
{
    auto it = m_patches.find(name);
    if (it == m_patches.end())
        throw std::out_of_range("UnstructuredBoundaryCondition: unknown patch '" + name + "'");
    return it->second;
}

bool UnstructuredBoundaryCondition::hasPatch(const std::string& name) const
{
    return m_patches.count(name) > 0;
}

void UnstructuredBoundaryCondition::applyVelocity(
    UnstructuredField<Eigen::Vector2d>& vel,
    const UnstructuredMesh&             mesh) const
{
    for (const auto& [name, patch] : m_patches)
    {
        if (patch.type != BoundaryType::INLET && patch.type != BoundaryType::WALL)
            continue;  // OUTLET/SYMMETRY: zero-gradient, leave unchanged

        // Patch may not exist on this mesh (e.g. config lists extra patches).
        if (mesh.numPatches() == 0)
            continue;

        bool found = false;
        for (int pi = 0; pi < mesh.numPatches(); ++pi)
        {
            if (mesh.patchName(pi) == name)
            {
                found = true;
                break;
            }
        }
        if (!found)
            continue;

        for (FaceId fid : mesh.patchFaces(name))
        {
            const MeshFace& f = mesh.face(fid);
            vel[f.owner] = patch.value;
        }
    }
}

void UnstructuredBoundaryCondition::applyPressure(
    UnstructuredField<double>&  p,
    const UnstructuredMesh&     mesh) const
{
    for (const auto& [name, patch] : m_patches)
    {
        if (patch.type != BoundaryType::OUTLET)
            continue;

        // Check patch exists in this mesh.
        bool found = false;
        for (int pi = 0; pi < mesh.numPatches(); ++pi)
        {
            if (mesh.patchName(pi) == name)
            {
                found = true;
                break;
            }
        }
        if (!found)
            continue;

        for (FaceId fid : mesh.patchFaces(name))
        {
            const MeshFace& f = mesh.face(fid);
            p[f.owner] = 0.0;
        }
    }
}

std::vector<CellId>
UnstructuredBoundaryCondition::collectCellsOfType(BoundaryType           type,
                                                   const UnstructuredMesh& mesh) const
{
    std::vector<CellId> cells;

    for (const auto& [name, patch] : m_patches)
    {
        if (patch.type != type)
            continue;

        bool found = false;
        for (int pi = 0; pi < mesh.numPatches(); ++pi)
        {
            if (mesh.patchName(pi) == name)
            {
                found = true;
                break;
            }
        }
        if (!found)
            continue;

        for (FaceId fid : mesh.patchFaces(name))
        {
            const MeshFace& f = mesh.face(fid);
            cells.push_back(f.owner);
        }
    }

    // Sort and deduplicate.
    std::sort(cells.begin(), cells.end(),
              [](CellId a, CellId b) { return toInt(a) < toInt(b); });
    cells.erase(std::unique(cells.begin(), cells.end()), cells.end());

    return cells;
}
