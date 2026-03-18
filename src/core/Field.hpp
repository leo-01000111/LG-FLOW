#pragma once

#include "core/Mesh.hpp"

#include <Eigen/Dense>
#include <stdexcept>
#include <vector>

/**
 * @brief Cell-centred scalar or vector field on a structured 2D mesh.
 *
 * Stores one value of type T per cell. Supported specialisations:
 *   - Field<double>          — scalar field (pressure, temperature, …)
 *   - Field<Eigen::Vector2d> — 2D vector field (velocity)
 *
 * @tparam T Value type stored at each cell (double or Eigen::Vector2d).
 */
template <typename T>
class Field
{
public:
    /**
     * @brief Constructs an uninitialised field associated with a mesh.
     * @param mesh Reference to the mesh that defines the topology.
     */
    explicit Field(const Mesh& mesh);

    /**
     * @brief Constructs a field and fills every cell with initialValue.
     * @param mesh         Reference to the mesh.
     * @param initialValue Value assigned to every cell at construction.
     */
    Field(const Mesh& mesh, const T& initialValue);

    /**
     * @brief Read-write access to the value at cell (i, j).
     * @param i x-direction index.
     * @param j y-direction index.
     * @return Reference to the cell value.
     */
    T& operator()(int i, int j);

    /**
     * @brief Read-only access to the value at cell (i, j).
     * @param i x-direction index.
     * @param j y-direction index.
     * @return Const reference to the cell value.
     */
    const T& operator()(int i, int j) const;

    /**
     * @brief Read-write access by linear cell index.
     * @param cellId Linear index = i * Ny + j.
     */
    T& operator[](int cellId);

    /** @brief Read-only access by linear cell index. */
    const T& operator[](int cellId) const;

    /**
     * @brief Sets every cell to the given value.
     * @param value The value to assign to all cells.
     */
    void setAll(const T& value);

    /**
     * @brief Computes the L2 norm of the field.
     *
     * For scalar fields: sqrt(Σ val²). For vector fields: sqrt(Σ |val|²).
     *
     * @return L2 norm (dimensionless ratio if field is normalised).
     */
    [[nodiscard]] double norm() const;

    /**
     * @brief Element-wise field addition.
     * @param other Field with the same mesh and type.
     * @return New field containing (*this)[i] + other[i] for each cell.
     */
    [[nodiscard]] Field<T> operator+(const Field<T>& other) const;

    /**
     * @brief Element-wise scalar multiplication.
     * @param scalar Multiplier applied to every cell value.
     * @return New field containing scalar * (*this)[i] for each cell.
     */
    [[nodiscard]] Field<T> operator*(double scalar) const;

    /** @brief Number of cells in the field (== mesh.numCells()). */
    [[nodiscard]] int size() const;

    /** @brief Reference to the mesh this field lives on. */
    [[nodiscard]] const Mesh& mesh() const;

private:
    const Mesh*  m_mesh;
    std::vector<T> m_data;
};

// ── Template implementation ───────────────────────────────────────────────────

template <typename T>
Field<T>::Field(const Mesh& mesh)
    : m_mesh(&mesh)
    , m_data(static_cast<std::size_t>(mesh.numCells()))
{
}

template <typename T>
Field<T>::Field(const Mesh& mesh, const T& initialValue)
    : m_mesh(&mesh)
    , m_data(static_cast<std::size_t>(mesh.numCells()), initialValue)
{
}

template <typename T>
T& Field<T>::operator()(int i, int j)
{
    return m_data[static_cast<std::size_t>(m_mesh->cellIndex(i, j))];
}

template <typename T>
const T& Field<T>::operator()(int i, int j) const
{
    return m_data[static_cast<std::size_t>(m_mesh->cellIndex(i, j))];
}

template <typename T>
T& Field<T>::operator[](int cellId)
{
    return m_data[static_cast<std::size_t>(cellId)];
}

template <typename T>
const T& Field<T>::operator[](int cellId) const
{
    return m_data[static_cast<std::size_t>(cellId)];
}

template <typename T>
void Field<T>::setAll(const T& value)
{
    for (auto& v : m_data)
        v = value;
}

template <typename T>
double Field<T>::norm() const
{
    // STUB: L2 norm will be implemented in Milestone 2
    return 0.0;
}

template <typename T>
Field<T> Field<T>::operator+(const Field<T>& other) const
{
    // STUB: element-wise addition will be implemented in Milestone 2
    return Field<T>(*m_mesh);
}

template <typename T>
Field<T> Field<T>::operator*(double /*scalar*/) const
{
    // STUB: scalar multiplication will be implemented in Milestone 2
    return Field<T>(*m_mesh);
}

template <typename T>
int Field<T>::size() const
{
    return static_cast<int>(m_data.size());
}

template <typename T>
const Mesh& Field<T>::mesh() const
{
    return *m_mesh;
}
