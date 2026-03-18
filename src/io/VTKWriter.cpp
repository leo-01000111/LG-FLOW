#include "io/VTKWriter.hpp"

#include <fstream>
#include <iomanip>
#include <stdexcept>

// Reference: VTK XML File Formats, version 0.1
//   https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf
//   Section 3 (XML Formats), UnstructuredGrid element type.
//
// Layout: one VTK_VERTEX (type=1) per mesh cell centre (z=0).
// Cell data: "pressure" scalar, "velocity" 3-component vector (Vx, Vy, 0).

void VTKWriter::write(const std::string&            filename,
                      const Mesh&                   mesh,
                      const Field<double>&          pressureField,
                      const Field<Eigen::Vector2d>& velocityField)
{
    // Validate mesh/field size compatibility before touching the filesystem.
    if (pressureField.size() != mesh.numCells())
        throw std::invalid_argument(
            "VTKWriter::write: pressure field size ("
            + std::to_string(pressureField.size())
            + ") does not match mesh cell count ("
            + std::to_string(mesh.numCells()) + ")");

    if (velocityField.size() != mesh.numCells())
        throw std::invalid_argument(
            "VTKWriter::write: velocity field size ("
            + std::to_string(velocityField.size())
            + ") does not match mesh cell count ("
            + std::to_string(mesh.numCells()) + ")");

    std::ofstream out(filename);
    if (!out.is_open())
        throw std::runtime_error(
            "VTKWriter::write: cannot create file: " + filename);

    const int Nc = mesh.numCells();
    const int Nx = mesh.Nx();
    const int Ny = mesh.Ny();

    out << std::scientific << std::setprecision(8);

    // ── XML header and opening tags ──────────────────────────────────────────
    out << "<?xml version=\"1.0\"?>\n";
    out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\""
           " byte_order=\"LittleEndian\">\n";
    out << "  <UnstructuredGrid>\n";
    out << "    <Piece NumberOfPoints=\"" << Nc
        << "\" NumberOfCells=\"" << Nc << "\">\n";

    // ── Points: one point per cell centre at z=0 ─────────────────────────────
    // Row-major traversal matches the cell indexing: cell(i,j) → i*Ny+j.
    out << "      <Points>\n";
    out << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\""
           " format=\"ascii\">\n";
    for (int i = 0; i < Nx; ++i)
    {
        for (int j = 0; j < Ny; ++j)
        {
            const Eigen::Vector2d c = mesh.getCellCenter(i, j);
            out << "          " << c.x() << " " << c.y() << " 0.0\n";
        }
    }
    out << "        </DataArray>\n";
    out << "      </Points>\n";

    // ── Cells: each cell is a VTK_VERTEX (type=1) referencing its own point ──
    out << "      <Cells>\n";

    // connectivity: point k belongs to cell k (identity mapping)
    out << "        <DataArray type=\"Int32\" Name=\"connectivity\""
           " format=\"ascii\">\n          ";
    for (int k = 0; k < Nc; ++k)
        out << " " << k;
    out << "\n        </DataArray>\n";

    // offsets: cell k ends at point index k+1
    out << "        <DataArray type=\"Int32\" Name=\"offsets\""
           " format=\"ascii\">\n          ";
    for (int k = 1; k <= Nc; ++k)
        out << " " << k;
    out << "\n        </DataArray>\n";

    // types: VTK_VERTEX = 1
    out << "        <DataArray type=\"UInt8\" Name=\"types\""
           " format=\"ascii\">\n          ";
    for (int k = 0; k < Nc; ++k)
        out << " 1";
    out << "\n        </DataArray>\n";

    out << "      </Cells>\n";

    // ── Cell data arrays ──────────────────────────────────────────────────────
    out << "      <CellData>\n";

    // Pressure (scalar [Pa])
    out << "        <DataArray type=\"Float64\" Name=\"pressure\""
           " NumberOfComponents=\"1\" format=\"ascii\">\n";
    for (int i = 0; i < Nx; ++i)
    {
        for (int j = 0; j < Ny; ++j)
            out << "          " << pressureField(i, j) << "\n";
    }
    out << "        </DataArray>\n";

    // Velocity (3-component vector: Vx, Vy, 0 [m/s])
    out << "        <DataArray type=\"Float64\" Name=\"velocity\""
           " NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (int i = 0; i < Nx; ++i)
    {
        for (int j = 0; j < Ny; ++j)
        {
            const Eigen::Vector2d& v = velocityField(i, j);
            out << "          " << v.x() << " " << v.y() << " 0.0\n";
        }
    }
    out << "        </DataArray>\n";

    out << "      </CellData>\n";
    out << "    </Piece>\n";
    out << "  </UnstructuredGrid>\n";
    out << "</VTKFile>\n";
}
