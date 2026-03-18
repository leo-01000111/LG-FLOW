#pragma once

#include "core/Field.hpp"
#include "core/Mesh.hpp"

#include <Eigen/Dense>
#include <string>

/**
 * @brief Writes simulation fields to VTK UnstructuredGrid (.vtu) format.
 *
 * Output is ASCII legacy VTK so it can be opened in ParaView without
 * additional libraries. Each call to write() produces one time-step snapshot.
 *
 * Field data is written as CELL_DATA arrays (cell-centred FVM convention).
 *
 * Reference: VTK File Formats — https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf
 */
class VTKWriter
{
public:
    VTKWriter() = default;

    /**
     * @brief Writes pressure and velocity fields to a .vtu file.
     *
     * Output fields:
     *   - "pressure"  — scalar, one value per cell [Pa]
     *   - "velocity"  — 3-component vector (Vx, Vy, 0) per cell [m/s]
     *
     * @param filename      Path to the output .vtu file (created or overwritten).
     * @param mesh          Mesh defining cell geometry.
     * @param pressureField Scalar pressure field.
     * @param velocityField Vector velocity field.
     * @throws std::runtime_error if the file cannot be created.
     */
    void write(const std::string&            filename,
               const Mesh&                   mesh,
               const Field<double>&          pressureField,
               const Field<Eigen::Vector2d>& velocityField);
};
