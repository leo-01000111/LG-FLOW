#include "core/Mesh.hpp"

#include <stdexcept>
#include <string>

// Uniform structured Cartesian mesh for 2D cell-centred FVM.
// Indexing spec (must match header contract):
//   cell(i,j)  → i*Ny + j,  i ∈ [0,Nx),  j ∈ [0,Ny)
//   x-face id  → i*Ny + j,  i ∈ [0,Nx],  j ∈ [0,Ny)   (xFaceCount = (Nx+1)*Ny)
//   y-face id  → xFaceCount + i*(Ny+1) + j,  i ∈ [0,Nx),  j ∈ [0,Ny]

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

    const int nCells     = Nx * Ny;
    const int xFaceCount = (Nx + 1) * Ny;
    const int nFaces     = xFaceCount + Nx * (Ny + 1);

    // ── Cell centres ─────────────────────────────────────────────────────────
    // Centre of cell(i,j) = ((i+0.5)*dx, (j+0.5)*dy)
    m_cellCentres.resize(static_cast<std::size_t>(nCells));
    for (int i = 0; i < Nx; ++i)
        for (int j = 0; j < Ny; ++j)
            m_cellCentres[static_cast<std::size_t>(i * Ny + j)] =
                {(i + 0.5) * m_dx, (j + 0.5) * m_dy};

    // ── Cell volumes (areas in 2-D) ───────────────────────────────────────────
    m_cellVolumes.assign(static_cast<std::size_t>(nCells), m_dx * m_dy);

    // ── Face geometry and connectivity ────────────────────────────────────────
    m_faceAreas.resize(static_cast<std::size_t>(nFaces));
    m_faceNormals.resize(static_cast<std::size_t>(nFaces));
    m_faceNeighbours.resize(static_cast<std::size_t>(nFaces));

    // X-faces: perpendicular to x, area = dy
    //   left  boundary (i=0)  : owner=cell(0,j),    neighbour=-1, normal=(-1, 0)
    //   right boundary (i=Nx) : owner=cell(Nx-1,j), neighbour=-1, normal=(+1, 0)
    //   interior (0<i<Nx)     : owner=cell(i-1,j),  neighbour=cell(i,j), normal=(+1, 0)
    for (int i = 0; i <= Nx; ++i)
    {
        for (int j = 0; j < Ny; ++j)
        {
            const auto fid = static_cast<std::size_t>(i * Ny + j);
            m_faceAreas[fid] = m_dy;

            if (i == 0)
            {
                m_faceNormals[fid]    = {-1.0, 0.0};
                m_faceNeighbours[fid] = {cellIndex(0, j), -1};
            }
            else if (i == Nx)
            {
                m_faceNormals[fid]    = {+1.0, 0.0};
                m_faceNeighbours[fid] = {cellIndex(Nx - 1, j), -1};
            }
            else
            {
                // Interior: normal points from owner (i-1) toward neighbour (i)
                m_faceNormals[fid]    = {+1.0, 0.0};
                m_faceNeighbours[fid] = {cellIndex(i - 1, j), cellIndex(i, j)};
            }
        }
    }

    // Y-faces: perpendicular to y, area = dx
    //   bottom boundary (j=0)  : owner=cell(i,0),    neighbour=-1, normal=(0,-1)
    //   top    boundary (j=Ny) : owner=cell(i,Ny-1), neighbour=-1, normal=(0,+1)
    //   interior (0<j<Ny)      : owner=cell(i,j-1),  neighbour=cell(i,j), normal=(0,+1)
    for (int i = 0; i < Nx; ++i)
    {
        for (int j = 0; j <= Ny; ++j)
        {
            const auto fid =
                static_cast<std::size_t>(xFaceCount + i * (Ny + 1) + j);
            m_faceAreas[fid] = m_dx;

            if (j == 0)
            {
                m_faceNormals[fid]    = {0.0, -1.0};
                m_faceNeighbours[fid] = {cellIndex(i, 0), -1};
            }
            else if (j == Ny)
            {
                m_faceNormals[fid]    = {0.0, +1.0};
                m_faceNeighbours[fid] = {cellIndex(i, Ny - 1), -1};
            }
            else
            {
                // Interior: normal points from owner (j-1) toward neighbour (j)
                m_faceNormals[fid]    = {0.0, +1.0};
                m_faceNeighbours[fid] = {cellIndex(i, j - 1), cellIndex(i, j)};
            }
        }
    }
}

Eigen::Vector2d Mesh::getCellCenter(int i, int j) const
{
    if (i < 0 || i >= m_Nx || j < 0 || j >= m_Ny)
        throw std::out_of_range(
            "Mesh::getCellCenter: index (" + std::to_string(i) + ", " +
            std::to_string(j) + ") out of range");
    return m_cellCentres[static_cast<std::size_t>(i * m_Ny + j)];
}

double Mesh::getFaceArea(int faceId) const
{
    if (faceId < 0 || faceId >= static_cast<int>(m_faceAreas.size()))
        throw std::out_of_range(
            "Mesh::getFaceArea: faceId " + std::to_string(faceId) + " out of range");
    return m_faceAreas[static_cast<std::size_t>(faceId)];
}

Eigen::Vector2d Mesh::getFaceNormal(int faceId) const
{
    if (faceId < 0 || faceId >= static_cast<int>(m_faceNormals.size()))
        throw std::out_of_range(
            "Mesh::getFaceNormal: faceId " + std::to_string(faceId) + " out of range");
    return m_faceNormals[static_cast<std::size_t>(faceId)];
}

std::array<int, 2> Mesh::getNeighbors(int faceId) const
{
    if (faceId < 0 || faceId >= static_cast<int>(m_faceNeighbours.size()))
        throw std::out_of_range(
            "Mesh::getNeighbors: faceId " + std::to_string(faceId) + " out of range");
    return m_faceNeighbours[static_cast<std::size_t>(faceId)];
}

double Mesh::getCellVolume(int cellId) const
{
    if (cellId < 0 || cellId >= static_cast<int>(m_cellVolumes.size()))
        throw std::out_of_range(
            "Mesh::getCellVolume: cellId " + std::to_string(cellId) + " out of range");
    return m_cellVolumes[static_cast<std::size_t>(cellId)];
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
