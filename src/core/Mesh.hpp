#pragma once

#include <Eigen/Dense>
#include <vector>
#include <array>

/**
 * @brief Structured 2D Cartesian mesh for cell-centered FVM.
 *
 * The mesh is an Nx × Ny grid of rectangular cells.
 * Cells are indexed in row-major order: cell(i, j) → linear index i*Ny + j,
 * where i ∈ [0, Nx) is the x-direction index and j ∈ [0, Ny) is the
 * y-direction index.
 *
 * Face convention (per cell):
 *   - face 0: west  (−x)
 *   - face 1: east  (+x)
 *   - face 2: south (−y)
 *   - face 3: north (+y)
 *
 * All geometric quantities are in SI units (metres).
 */
class Mesh
{
public:
    Mesh() = default;

    /**
     * @brief Builds the mesh topology and geometry.
     *
     * Uniform Cartesian grid spanning [0, Lx] × [0, Ly] with Nx × Ny cells.
     *
     * @param Nx Number of cells in the x-direction.
     * @param Ny Number of cells in the y-direction.
     * @param Lx Domain length in the x-direction [m].
     * @param Ly Domain length in the y-direction [m].
     * @throws std::invalid_argument if Nx or Ny is zero, or Lx/Ly is non-positive.
     */
    void load(int Nx, int Ny, double Lx, double Ly);

    /**
     * @brief Returns the physical coordinates of the centre of cell (i, j).
     * @param i  x-direction cell index ∈ [0, Nx).
     * @param j  y-direction cell index ∈ [0, Ny).
     * @return 2D position vector [m].
     */
    [[nodiscard]] Eigen::Vector2d getCellCenter(int i, int j) const;

    /**
     * @brief Returns the face area (length in 2D) for the given global face ID.
     *
     * Face IDs are assigned in x-direction sweeps first, then y-direction.
     *
     * @param faceId Global face index.
     * @return Face length [m].
     */
    [[nodiscard]] double getFaceArea(int faceId) const;

    /**
     * @brief Returns the outward-pointing unit normal for the given face ID.
     * @param faceId Global face index.
     * @return Unit normal vector (dimensionless).
     */
    [[nodiscard]] Eigen::Vector2d getFaceNormal(int faceId) const;

    /**
     * @brief Returns the two cell indices that share the given interior face.
     *
     * For boundary faces, the second index is -1 (no neighbour).
     *
     * @param faceId Global face index.
     * @return Array of two cell indices {owner, neighbour}.
     */
    [[nodiscard]] std::array<int, 2> getNeighbors(int faceId) const;

    /**
     * @brief Returns the volume (area in 2D) of a cell.
     * @param cellId Linear cell index.
     * @return Cell volume [m²].
     */
    [[nodiscard]] double getCellVolume(int cellId) const;

    /** @brief Total number of cells (Nx * Ny). */
    [[nodiscard]] int numCells() const;

    /** @brief Total number of faces (interior + boundary). */
    [[nodiscard]] int numFaces() const;

    /** @brief Number of cells in the x-direction. */
    [[nodiscard]] int Nx() const;

    /** @brief Number of cells in the y-direction. */
    [[nodiscard]] int Ny() const;

    /**
     * @brief Converts a 2D index pair to a linear cell index.
     * @param i x-direction index.
     * @param j y-direction index.
     * @return Linear index = i * m_Ny + j.
     */
    [[nodiscard]] int cellIndex(int i, int j) const;

private:
    int    m_Nx{0};
    int    m_Ny{0};
    double m_dx{0.0};  ///< Uniform cell width  [m]
    double m_dy{0.0};  ///< Uniform cell height [m]

    /// Cell-centre coordinates, size m_Nx * m_Ny
    std::vector<Eigen::Vector2d> m_cellCentres;

    /// Cell volumes (areas in 2D), size m_Nx * m_Ny
    std::vector<double> m_cellVolumes;

    /// Face areas, size numFaces
    std::vector<double> m_faceAreas;

    /// Outward unit normals per face
    std::vector<Eigen::Vector2d> m_faceNormals;

    /// Owner-neighbour connectivity per face
    std::vector<std::array<int, 2>> m_faceNeighbours;
};
