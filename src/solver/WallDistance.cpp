#include "solver/WallDistance.hpp"

#include <limits>
#include <queue>
#include <vector>

Field<double> WallDistance::compute(const Mesh& mesh, const BoundaryCondition& bc)
{
    const int    N  = mesh.numCells();
    const int    Nx = mesh.Nx();
    const int    Ny = mesh.Ny();

    // Seed distance for wall-adjacent cells: half the minimum cell dimension
    // equals the distance from the cell centre to the nearest wall face.
    const double dx      = mesh.getCellCenter(1, 0).x() - mesh.getCellCenter(0, 0).x();
    const double dy      = mesh.getCellCenter(0, 1).y() - mesh.getCellCenter(0, 0).y();
    const double halfMin = 0.5 * std::min(dx, dy);

    Field<double>      wallDist(mesh, std::numeric_limits<double>::max());
    std::vector<bool>  inQueue(static_cast<std::size_t>(N), false);
    std::queue<int>    bfsQueue;

    // Seed all WALL cells with the half-cell distance and add them to the queue.
    for (int c : bc.collectCellsOfType(BoundaryType::WALL, mesh))
    {
        wallDist[c] = halfMin;
        if (!inQueue[static_cast<std::size_t>(c)])
        {
            inQueue[static_cast<std::size_t>(c)] = true;
            bfsQueue.push(c);
        }
    }

    // BFS: propagate distances to interior cells.
    // For each cell popped from the queue, update its four structured neighbors
    // (i±1, j±1) with the candidate distance: dist[c] + |centre[n] - centre[c]|.
    // If the candidate is smaller than the current stored distance, update and re-queue.
    while (!bfsQueue.empty())
    {
        const int c = bfsQueue.front();
        bfsQueue.pop();
        inQueue[static_cast<std::size_t>(c)] = false;

        const int io = c / Ny;
        const int jo = c % Ny;
        const Eigen::Vector2d xo = mesh.getCellCenter(io, jo);

        // Neighbor offsets: west, east, south, north
        const int di[4] = {-1,  1,  0,  0};
        const int dj[4] = { 0,  0, -1,  1};

        for (int dir = 0; dir < 4; ++dir)
        {
            const int ni = io + di[dir];
            const int nj = jo + dj[dir];

            if (ni < 0 || ni >= Nx || nj < 0 || nj >= Ny)
                continue;

            const int n = ni * Ny + nj;
            const Eigen::Vector2d xn = mesh.getCellCenter(ni, nj);
            const double candidate = wallDist[c] + (xn - xo).norm();

            if (candidate < wallDist[n])
            {
                wallDist[n] = candidate;
                // Only push if not already in the queue to avoid redundant work.
                if (!inQueue[static_cast<std::size_t>(n)])
                {
                    inQueue[static_cast<std::size_t>(n)] = true;
                    bfsQueue.push(n);
                }
            }
        }
    }

    return wallDist;
}
