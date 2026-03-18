#include "core/Mesh.hpp"

#include <stdexcept>

// STUB: Full mesh construction (cell centres, face connectivity, normals) will
// be implemented in Milestone 2.

void Mesh::load(int Nx, int Ny, double Lx, double Ly)
{
    if (Nx <= 0 || Ny <= 0)
        throw std::invalid_argument("Mesh::load: Nx and Ny must be positive");
    if (Lx <= 0.0 || Ly <= 0.0)
        throw std::invalid_argument("Mesh::load: Lx and Ly must be positive");

    m_Nx = Nx;
    m_Ny = Ny;
    m_dx = Lx / static_cast<double>(Nx);
    m_dy = Ly / static_cast<double>(Ny);

    // STUB: populate m_cellCentres, m_cellVolumes, m_faceAreas,
    // m_faceNormals, m_faceNeighbours
}

Eigen::Vector2d Mesh::getCellCenter(int /*i*/, int /*j*/) const
{
    // STUB: returns cell centre at (i + 0.5)*dx, (j + 0.5)*dy
    return Eigen::Vector2d::Zero();
}

double Mesh::getFaceArea(int /*faceId*/) const
{
    // STUB: returns m_dx or m_dy depending on face orientation
    return 0.0;
}

Eigen::Vector2d Mesh::getFaceNormal(int /*faceId*/) const
{
    // STUB: returns outward unit normal for the requested face
    return Eigen::Vector2d::Zero();
}

std::array<int, 2> Mesh::getNeighbors(int /*faceId*/) const
{
    // STUB: returns {owner, neighbour}; neighbour = -1 on boundary
    return {-1, -1};
}

double Mesh::getCellVolume(int /*cellId*/) const
{
    // STUB: returns m_dx * m_dy for uniform grid
    return 0.0;
}

int Mesh::numCells() const { return m_Nx * m_Ny; }

int Mesh::numFaces() const
{
    // x-faces: (Nx+1)*Ny, y-faces: Nx*(Ny+1)
    return (m_Nx + 1) * m_Ny + m_Nx * (m_Ny + 1);
}

int Mesh::Nx() const { return m_Nx; }
int Mesh::Ny() const { return m_Ny; }

int Mesh::cellIndex(int i, int j) const { return i * m_Ny + j; }
