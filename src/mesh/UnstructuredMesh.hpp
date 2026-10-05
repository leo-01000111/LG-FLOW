#pragma once

#include <Eigen/Dense>
#include <array>
#include <span>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

// ── Strong-typed ID wrappers ──────────────────────────────────────────────────
// Compile-time protection against mixing node/face/cell indices.

/** @brief Strong-typed ID for mesh nodes. Underlying type int. */
enum class NodeId : int {};

/** @brief Strong-typed ID for mesh faces. Underlying type int. */
enum class FaceId : int {};

/** @brief Strong-typed ID for mesh cells. Underlying type int. */
enum class CellId : int {};

/** @brief Converts NodeId to int for vector indexing. */
inline int toInt(NodeId id) { return static_cast<int>(id); }

/** @brief Converts FaceId to int for vector indexing. */
inline int toInt(FaceId id) { return static_cast<int>(id); }

/** @brief Converts CellId to int for vector indexing. */
inline int toInt(CellId id) { return static_cast<int>(id); }

// ── Hash for std::pair<int,int> ────────────────────────────────────────────────
/**
 * @brief Hash functor for std::pair<int,int> used in unordered_map edge keys.
 *
 * Combines the two ints into a single 64-bit key using a large prime multiplier.
 * Assumes node IDs are well below 1 000 000.
 */
struct PairHash
{
    std::size_t operator()(const std::pair<int, int>& p) const
    {
        return std::hash<long long>{}(static_cast<long long>(p.first) * 1000000LL + p.second);
    }
};

// ── Plain data types ──────────────────────────────────────────────────────────

/** @brief A mesh node with a 2D position. */
struct MeshNode
{
    Eigen::Vector2d position; ///< Node coordinates [m]
};

/**
 * @brief A mesh face (edge in 2D) connecting two nodes and separating two cells.
 *
 * Owner: the cell whose outward normal convention the stored normal follows.
 * Neighbour: CellId{-1} for boundary faces.
 */
struct MeshFace
{
    std::array<NodeId, 2> nodes;      ///< Endpoint node IDs
    CellId                owner;      ///< Owning cell (normal points away from owner interior)
    CellId                neighbour;  ///< Adjacent cell; CellId{-1} for boundary faces
    Eigen::Vector2d       normal;     ///< Outward unit normal from owner's perspective
    double                area;       ///< Face length [m]
    Eigen::Vector2d       centre;     ///< Face midpoint [m]
    int                   patchId;    ///< Index into patch list; -1 for interior faces
};

/**
 * @brief A mesh cell (polygon) defined by an ordered CCW node list.
 *
 * Faces list is populated by buildConnectivity().
 */
struct MeshCell
{
    std::vector<NodeId>  nodes;   ///< Node IDs in CCW order
    std::vector<FaceId>  faces;   ///< Bounding face IDs
    Eigen::Vector2d      centre;  ///< Cell centroid [m]
    double               volume;  ///< Cell area [m²]
};

// ── UnstructuredMesh ──────────────────────────────────────────────────────────

/**
 * @brief Unstructured 2D cell-centred FVM mesh.
 *
 * Supports triangular and quadrilateral cells. Face connectivity is built
 * from cell node lists by buildConnectivity(). Interior faces are stored
 * before boundary faces, grouped by patch.
 *
 * Lifecycle:
 *   1. Add nodes via addNode().
 *   2. Add cells via addCell().
 *   3. Register patches via addPatch().
 *   4. Call buildConnectivity() with the boundary-edge-to-patch map.
 */
class UnstructuredMesh
{
public:
    UnstructuredMesh() = default;

    /** @brief Total number of nodes. */
    [[nodiscard]] int numNodes() const;

    /** @brief Total number of faces (interior + boundary). */
    [[nodiscard]] int numFaces() const;

    /** @brief Total number of cells. */
    [[nodiscard]] int numCells() const;

    /** @brief Number of interior (non-boundary) faces. */
    [[nodiscard]] int numInteriorFaces() const;

    /**
     * @brief Returns the node with the given ID.
     * @param id Node ID.
     */
    [[nodiscard]] const MeshNode& node(NodeId id) const;

    /**
     * @brief Returns the face with the given ID.
     * @param id Face ID.
     */
    [[nodiscard]] const MeshFace& face(FaceId id) const;

    /**
     * @brief Returns the cell with the given ID.
     * @param id Cell ID.
     */
    [[nodiscard]] const MeshCell& cell(CellId id) const;

    /** @brief Total number of registered boundary patches. */
    [[nodiscard]] int numPatches() const;

    /**
     * @brief Returns the name of patch at the given index.
     * @param patchIndex Patch index in [0, numPatches()).
     */
    [[nodiscard]] const std::string& patchName(int patchIndex) const;

    /**
     * @brief Returns the face IDs belonging to a named patch.
     * @param name Patch name.
     * @throws std::out_of_range if name is not registered.
     */
    [[nodiscard]] const std::vector<FaceId>& patchFaces(const std::string& name) const;

    /**
     * @brief Returns the face IDs belonging to the patch at the given index.
     * @param patchIndex Patch index.
     */
    [[nodiscard]] const std::vector<FaceId>& patchFaces(int patchIndex) const;

    /**
     * @brief Returns true if the face is on a boundary (no neighbour cell).
     * @param f Face ID.
     */
    [[nodiscard]] bool isBoundaryFace(FaceId f) const;

    // ── Mutable interface (called only by GmshReader) ─────────────────────────

    /**
     * @brief Adds a node and returns its ID.
     * @param n Node to add.
     */
    NodeId addNode(const MeshNode& n);

    /**
     * @brief Adds a cell defined by an ordered CCW node ID list.
     *
     * Centroid and volume are computed immediately from the node positions.
     *
     * @param nodeIds CCW-ordered node IDs of the cell polygon.
     * @return The new cell's ID.
     */
    CellId addCell(std::span<const NodeId> nodeIds);

    /**
     * @brief Registers a patch name and returns its index.
     *
     * If the name is already registered the existing index is returned.
     *
     * @param name Patch name (e.g. "inlet", "wall_bottom").
     * @return Patch index.
     */
    int addPatch(const std::string& name);

    /**
     * @brief Builds face connectivity from cell node lists.
     *
     * Algorithm:
     *   1. For each cell, iterate edges (consecutive node pairs, wrapping).
     *   2. Canonical edge key: {min(n0,n1), max(n0,n1)} — shared face detection.
     *   3. First encounter: create face with owner = cell. Second: assign neighbour.
     *   4. Faces with no neighbour are boundary faces.
     *   5. Interior faces first, then boundary faces sorted by patchId.
     *   6. Rebuild each cell's face list from the final ordering.
     *
     * @param boundaryEdgeToPatch Maps canonical edge key to patch index.
     * @throws std::runtime_error if an edge is shared by 3+ cells.
     */
    void buildConnectivity(
        const std::unordered_map<std::pair<int, int>, int, PairHash>& boundaryEdgeToPatch);

private:
    std::vector<MeshNode>  m_nodes;
    std::vector<MeshFace>  m_faces;
    std::vector<MeshCell>  m_cells;
    std::vector<std::string>         m_patchNames;
    std::vector<std::vector<FaceId>> m_patchFaces;
    std::unordered_map<std::string, int> m_patchIndex;
    int m_numInteriorFaces{0};

    /**
     * @brief Computes the centroid and area of a polygon from its node list.
     *
     * Uses the shoelace formula for a CCW polygon.
     * Reference: computational geometry; Gauss area formula.
     *
     * @param nodeIds Ordered node IDs.
     * @param nodes   Node position array.
     * @return {centroid, area}
     */
    [[nodiscard]] static std::pair<Eigen::Vector2d, double>
    computeCentroidAndArea(std::span<const NodeId>      nodeIds,
                           const std::vector<MeshNode>& nodes);
};
