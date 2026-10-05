#pragma once

#include "mesh/UnstructuredMesh.hpp"

#include <Eigen/Dense>
#include <cmath>
#include <stdexcept>
#include <type_traits>
#include <vector>

/**
 * @brief Cell-centred scalar or vector field on an unstructured 2D mesh.
 *
 * Stores one value of type T per cell, indexed by CellId.
 * Mirrors the interface of Field<T> (structured), but uses unstructured
 * mesh topology and CellId-based access.
 *
 * @tparam T Value type (double or Eigen::Vector2d).
 */
template <typename T>
class UnstructuredField
{
public:
    /**
     * @brief Constructs an uninitialised field associated with a mesh.
     * @param mesh Reference to the unstructured mesh (must outlive this object).
     */
    explicit UnstructuredField(const UnstructuredMesh& mesh)
        : m_mesh(&mesh)
        , m_data(static_cast<std::size_t>(mesh.numCells()))
    {
    }

    /**
     * @brief Constructs a field and fills every cell with initialValue.
     * @param mesh         Reference to the mesh.
     * @param initialValue Value assigned to every cell at construction.
     */
    UnstructuredField(const UnstructuredMesh& mesh, const T& initialValue)
        : m_mesh(&mesh)
        , m_data(static_cast<std::size_t>(mesh.numCells()), initialValue)
    {
    }

    /**
     * @brief Read-write access to the value at the given cell.
     * @param id Cell ID.
     */
    T& operator[](CellId id)
    {
        return m_data[static_cast<std::size_t>(toInt(id))];
    }

    /**
     * @brief Read-only access to the value at the given cell.
     * @param id Cell ID.
     */
    [[nodiscard]] const T& operator[](CellId id) const
    {
        return m_data[static_cast<std::size_t>(toInt(id))];
    }

    /**
     * @brief Sets every cell to the given value.
     * @param value Value to assign.
     */
    void setAll(const T& value)
    {
        for (auto& v : m_data)
            v = value;
    }

    /**
     * @brief Computes the L2 norm of the field.
     *
     * Scalar: sqrt(Σ val²). Vector: sqrt(Σ |val|²).
     *
     * @return L2 norm.
     */
    [[nodiscard]] double norm() const
    {
        double sum = 0.0;
        for (const auto& v : m_data)
        {
            if constexpr (std::is_same_v<T, double>)
                sum += v * v;
            else
                sum += v.squaredNorm();
        }
        return std::sqrt(sum);
    }

    /** @brief Number of cells (== mesh.numCells()). */
    [[nodiscard]] int size() const
    {
        return static_cast<int>(m_data.size());
    }

    /** @brief Reference to the mesh this field lives on. */
    [[nodiscard]] const UnstructuredMesh& mesh() const
    {
        return *m_mesh;
    }

    /**
     * @brief Element-wise field addition.
     * @param other Field on the same mesh.
     * @return New field with (*this)[c] + other[c].
     * @throws std::invalid_argument if mesh pointers differ.
     */
    [[nodiscard]] UnstructuredField<T> operator+(const UnstructuredField<T>& other) const
    {
        if (m_mesh != other.m_mesh || size() != other.size())
            throw std::invalid_argument(
                "UnstructuredField::operator+: fields must be on the same mesh");

        UnstructuredField<T> result(*m_mesh);
        const std::size_t n = m_data.size();
        for (std::size_t k = 0; k < n; ++k)
            result.m_data[k] = m_data[k] + other.m_data[k];
        return result;
    }

    /**
     * @brief Element-wise scalar multiplication.
     * @param scalar Multiplier applied to every cell value.
     * @return New field with scalar * (*this)[c].
     */
    [[nodiscard]] UnstructuredField<T> operator*(double scalar) const
    {
        UnstructuredField<T> result(*m_mesh);
        const std::size_t n = m_data.size();
        for (std::size_t k = 0; k < n; ++k)
            result.m_data[k] = m_data[k] * scalar;
        return result;
    }

private:
    const UnstructuredMesh* m_mesh;
    std::vector<T>          m_data;
};
