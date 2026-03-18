#pragma once

#include "core/Mesh.hpp"

#include <string>

/**
 * @brief Reads a structured mesh from a simple text-based format.
 *
 * File format (subject to change in Milestone 2):
 * @code
 *   # LG-Flow structured mesh
 *   Nx Ny
 *   Lx Ly
 * @endcode
 *
 * Physical dimensions are in metres.
 */
class MeshReader
{
public:
    MeshReader() = default;

    /**
     * @brief Reads mesh dimensions from a file and populates a Mesh object.
     *
     * @param filepath Path to the mesh description file.
     * @param mesh     Mesh object to populate (existing data will be overwritten).
     * @throws std::runtime_error if the file cannot be opened or is malformed.
     */
    [[nodiscard]] bool read(const std::string& filepath, Mesh& mesh);
};
