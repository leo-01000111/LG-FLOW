#include "mesh/UnstructuredMesh.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>

// ── Accessors ──────────────────────────────────────────────────────────────────

int UnstructuredMesh::numNodes() const
{
    return static_cast<int>(m_nodes.size());
}

int UnstructuredMesh::numFaces() const
{
    return static_cast<int>(m_faces.size());
}

int UnstructuredMesh::numCells() const
{
    return static_cast<int>(m_cells.size());
}

int UnstructuredMesh::numInteriorFaces() const
{
    return m_numInteriorFaces;
}

const MeshNode& UnstructuredMesh::node(NodeId id) const
{
    return m_nodes[static_cast<std::size_t>(toInt(id))];
}

const MeshFace& UnstructuredMesh::face(FaceId id) const
{
    return m_faces[static_cast<std::size_t>(toInt(id))];
}

const MeshCell& UnstructuredMesh::cell(CellId id) const
{
    return m_cells[static_cast<std::size_t>(toInt(id))];
}

int UnstructuredMesh::numPatches() const
{
    return static_cast<int>(m_patchNames.size());
}

const std::string& UnstructuredMesh::patchName(int patchIndex) const
{
    return m_patchNames[static_cast<std::size_t>(patchIndex)];
}

const std::vector<FaceId>& UnstructuredMesh::patchFaces(const std::string& name) const
{
    auto it = m_patchIndex.find(name);
    if (it == m_patchIndex.end())
        throw std::out_of_range("UnstructuredMesh: unknown patch '" + name + "'");
    return m_patchFaces[static_cast<std::size_t>(it->second)];
}

const std::vector<FaceId>& UnstructuredMesh::patchFaces(int patchIndex) const
{
    return m_patchFaces[static_cast<std::size_t>(patchIndex)];
}

bool UnstructuredMesh::isBoundaryFace(FaceId f) const
{
    return toInt(m_faces[static_cast<std::size_t>(toInt(f))].neighbour) < 0;
}

// ── Mutable interface ─────────────────────────────────────────────────────────

NodeId UnstructuredMesh::addNode(const MeshNode& n)
{
    const auto id = static_cast<int>(m_nodes.size());
    m_nodes.push_back(n);
    return NodeId{id};
}

CellId UnstructuredMesh::addCell(std::span<const NodeId> nodeIds)
{
    const auto id = static_cast<int>(m_cells.size());

    MeshCell c;
    c.nodes.assign(nodeIds.begin(), nodeIds.end());

    auto [ctr, vol] = computeCentroidAndArea(nodeIds, m_nodes);
    c.centre = ctr;
    c.volume = vol;

    m_cells.push_back(std::move(c));
    return CellId{id};
}

int UnstructuredMesh::addPatch(const std::string& name)
{
    auto it = m_patchIndex.find(name);
    if (it != m_patchIndex.end())
        return it->second;

    const int idx = static_cast<int>(m_patchNames.size());
    m_patchNames.push_back(name);
    m_patchFaces.emplace_back();
    m_patchIndex[name] = idx;
    return idx;
}

// ── Centroid and area ─────────────────────────────────────────────────────────

std::pair<Eigen::Vector2d, double>
UnstructuredMesh::computeCentroidAndArea(std::span<const NodeId>      nodeIds,
                                         const std::vector<MeshNode>& nodes)
{
    // Shoelace formula for CCW polygon area and centroid.
    // Reference: Computational Geometry — Gauss area formula.
    const int n = static_cast<int>(nodeIds.size());
    double area = 0.0;
    Eigen::Vector2d centroid = Eigen::Vector2d::Zero();

    for (int i = 0; i < n; ++i)
    {
        const int j = (i + 1) % n;
        const Eigen::Vector2d& pi = nodes[static_cast<std::size_t>(toInt(nodeIds[i]))].position;
        const Eigen::Vector2d& pj = nodes[static_cast<std::size_t>(toInt(nodeIds[j]))].position;

        const double cross = pi.x() * pj.y() - pj.x() * pi.y();
        area     += cross;
        centroid += (pi + pj) * cross;
    }

    area     *= 0.5;
    if (std::abs(area) < 1e-30)
    {
        // Degenerate cell — return geometric average as fallback.
        Eigen::Vector2d avg = Eigen::Vector2d::Zero();
        for (const auto& id : nodeIds)
            avg += nodes[static_cast<std::size_t>(toInt(id))].position;
        avg /= static_cast<double>(n);
        return {avg, 0.0};
    }
    centroid /= (6.0 * area);  // complete shoelace centroid formula
    return {centroid, std::abs(area)};
}

// ── buildConnectivity ─────────────────────────────────────────────────────────

void UnstructuredMesh::buildConnectivity(
    const std::unordered_map<std::pair<int, int>, int, PairHash>& boundaryEdgeToPatch)
{
    // Phase 1: walk every cell's edge list and collect half-faces.
    // A "half-face" is one directed edge from a single cell's perspective.
    // When two cells share an edge, we combine them into one interior face.

    // edgeToFaceIdx: canonical edge key → index into rawFaces (temporary storage).
    std::unordered_map<std::pair<int, int>, int, PairHash> edgeToFaceIdx;

    // Temporary face storage before final sort.
    std::vector<MeshFace> rawFaces;
    rawFaces.reserve(static_cast<std::size_t>(m_cells.size()) * 3);

    // For each raw face, track which cells reference it (for face-cell rebuild).
    // We rebuild cell.faces at the end after final face ordering.
    // ownerFor[rawFaceIdx] and neighbourFor[rawFaceIdx] hold the cell IDs.

    const int numCells = static_cast<int>(m_cells.size());
    for (int ci = 0; ci < numCells; ++ci)
    {
        const MeshCell& c = m_cells[static_cast<std::size_t>(ci)];
        const int numNodes = static_cast<int>(c.nodes.size());

        for (int ei = 0; ei < numNodes; ++ei)
        {
            const int n0 = toInt(c.nodes[static_cast<std::size_t>(ei)]);
            const int n1 = toInt(c.nodes[static_cast<std::size_t>((ei + 1) % numNodes)]);

            // Canonical key: smaller node ID first.
            const std::pair<int, int> key = (n0 < n1)
                ? std::pair<int, int>{n0, n1}
                : std::pair<int, int>{n1, n0};

            auto [it, inserted] = edgeToFaceIdx.emplace(key, static_cast<int>(rawFaces.size()));

            if (inserted)
            {
                // First encounter — create face with this cell as owner.
                MeshFace f{};
                f.nodes[0]   = NodeId{n0};
                f.nodes[1]   = NodeId{n1};
                f.owner      = CellId{ci};
                f.neighbour  = CellId{-1};
                f.patchId    = -1;

                // Look up patch for this edge.
                auto pit = boundaryEdgeToPatch.find(key);
                if (pit != boundaryEdgeToPatch.end())
                    f.patchId = pit->second;

                rawFaces.push_back(f);
            }
            else
            {
                // Second encounter — assign this cell as neighbour.
                const int fIdx = it->second;
                if (toInt(rawFaces[static_cast<std::size_t>(fIdx)].neighbour) >= 0)
                    throw std::runtime_error(
                        "UnstructuredMesh::buildConnectivity: edge shared by 3+ cells "
                        "(nodes " + std::to_string(key.first) + ", " + std::to_string(key.second) + ")");

                rawFaces[static_cast<std::size_t>(fIdx)].neighbour = CellId{ci};
                // Interior faces have patchId -1 regardless of boundary map.
                rawFaces[static_cast<std::size_t>(fIdx)].patchId = -1;
            }
        }
    }

    // Phase 2: separate interior and boundary faces, then sort boundary by patchId.
    std::vector<MeshFace> interiorFaces;
    std::vector<MeshFace> boundaryFaces;
    interiorFaces.reserve(rawFaces.size());
    boundaryFaces.reserve(rawFaces.size());

    for (auto& f : rawFaces)
    {
        if (toInt(f.neighbour) < 0)
            boundaryFaces.push_back(std::move(f));
        else
            interiorFaces.push_back(std::move(f));
    }

    // Sort by index rather than by value to avoid aligned_storage issues with
    // Eigen::Vector2d members under MSVC's stable_sort implementation.
    {
        std::vector<std::size_t> order(boundaryFaces.size());
        for (std::size_t k = 0; k < order.size(); ++k) order[k] = k;
        std::stable_sort(order.begin(), order.end(),
            [&](std::size_t a, std::size_t b) {
                return boundaryFaces[a].patchId < boundaryFaces[b].patchId;
            });
        std::vector<MeshFace> sorted;
        sorted.reserve(boundaryFaces.size());
        for (std::size_t k : order)
            sorted.push_back(std::move(boundaryFaces[k]));
        boundaryFaces = std::move(sorted);
    }

    m_numInteriorFaces = static_cast<int>(interiorFaces.size());

    // Phase 3: merge into m_faces (interior first, then boundary).
    m_faces.clear();
    m_faces.reserve(interiorFaces.size() + boundaryFaces.size());
    for (auto& f : interiorFaces)
        m_faces.push_back(std::move(f));
    for (auto& f : boundaryFaces)
        m_faces.push_back(std::move(f));

    // Phase 4: compute face geometry (normal, area, centre) for each face.
    // Normal convention: outward from owner cell.
    // For edge (A→B): right-hand perpendicular is (dy, -dx) in 2D.
    // We pick the perpendicular that points away from the owner centroid.
    for (auto& f : m_faces)
    {
        const Eigen::Vector2d& posA = m_nodes[static_cast<std::size_t>(toInt(f.nodes[0]))].position;
        const Eigen::Vector2d& posB = m_nodes[static_cast<std::size_t>(toInt(f.nodes[1]))].position;

        const Eigen::Vector2d edge = posB - posA;
        f.area   = edge.norm();
        f.centre = 0.5 * (posA + posB);

        // Two candidate normals — choose the one pointing away from owner centroid.
        // Right-hand perp of (dx, dy) is (dy, -dx); left-hand is (-dy, dx).
        Eigen::Vector2d n{edge.y(), -edge.x()};

        if (f.area > 1e-30)
            n /= f.area;  // normalise

        const Eigen::Vector2d ownerCentre = m_cells[static_cast<std::size_t>(toInt(f.owner))].centre;
        // If normal points toward owner centroid, flip it.
        if (n.dot(f.centre - ownerCentre) < 0.0)
            n = -n;

        f.normal = n;
    }

    // Phase 5: rebuild each cell's face list from the final face ordering.
    // Build a lookup: canonical edge key → final FaceId.
    std::unordered_map<std::pair<int, int>, FaceId, PairHash> edgeToFinalFace;
    edgeToFinalFace.reserve(m_faces.size());

    for (int fi = 0; fi < static_cast<int>(m_faces.size()); ++fi)
    {
        const MeshFace& f = m_faces[static_cast<std::size_t>(fi)];
        const int n0 = toInt(f.nodes[0]);
        const int n1 = toInt(f.nodes[1]);
        const std::pair<int, int> key = (n0 < n1)
            ? std::pair<int, int>{n0, n1}
            : std::pair<int, int>{n1, n0};
        edgeToFinalFace.emplace(key, FaceId{fi});
    }

    for (auto& c : m_cells)
    {
        c.faces.clear();
        const int numCellNodes = static_cast<int>(c.nodes.size());
        c.faces.reserve(static_cast<std::size_t>(numCellNodes));

        for (int ei = 0; ei < numCellNodes; ++ei)
        {
            const int n0 = toInt(c.nodes[static_cast<std::size_t>(ei)]);
            const int n1 = toInt(c.nodes[static_cast<std::size_t>((ei + 1) % numCellNodes)]);

            const std::pair<int, int> key = (n0 < n1)
                ? std::pair<int, int>{n0, n1}
                : std::pair<int, int>{n1, n0};

            c.faces.push_back(edgeToFinalFace.at(key));
        }
    }

    // Phase 6: populate patch face lists.
    m_patchFaces.assign(m_patchNames.size(), std::vector<FaceId>{});

    for (int fi = 0; fi < static_cast<int>(m_faces.size()); ++fi)
    {
        const MeshFace& f = m_faces[static_cast<std::size_t>(fi)];
        if (f.patchId >= 0 && static_cast<std::size_t>(f.patchId) < m_patchFaces.size())
            m_patchFaces[static_cast<std::size_t>(f.patchId)].push_back(FaceId{fi});
    }
}
