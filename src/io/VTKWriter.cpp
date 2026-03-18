#include "io/VTKWriter.hpp"

void VTKWriter::write(const std::string&            /*filename*/,
                      const Mesh&                   /*mesh*/,
                      const Field<double>&          /*pressureField*/,
                      const Field<Eigen::Vector2d>& /*velocityField*/)
{
    // STUB: ASCII .vtu output will be implemented in Milestone 2.
    //
    // Implementation plan:
    //   1. Write VTK header and DATASET UNSTRUCTURED_GRID
    //   2. Write POINTS (cell-centre x, y, 0 for each cell)
    //   3. Write CELLS connectivity (each cell = 1 point, type VTK_VERTEX)
    //   4. Write CELL_DATA section:
    //      - SCALARS pressure float 1
    //      - VECTORS velocity float
}
