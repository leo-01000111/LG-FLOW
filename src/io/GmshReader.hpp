#pragma once

#include "mesh/UnstructuredMesh.hpp"

#include <istream>
#include <string>
#include <unordered_map>
#include <utility>

/**
 * @brief Reads Gmsh .msh format 2 (ASCII) files into an UnstructuredMesh.
 *
 * Supported element types:
 *   - Type 1: 2-node line (boundary edge → patch)
 *   - Type 2: 3-node triangle (fluid cell)
 *   - Type 3: 4-node quad (fluid cell)
 *
 * Physical groups:
 *   - dim=1 → boundary patch (name maps to patch in mesh)
 *   - dim=2 → fluid region (ignored; only one region for 2D)
 *
 * Usage:
 * @code
 *   UnstructuredMesh mesh;
 *   GmshReader reader;
 *   reader.read("channel.msh", mesh);
 * @endcode
 */
class GmshReader
{
public:
    GmshReader() = default;

    /**
     * @brief Reads a Gmsh .msh v2 file and populates the given mesh.
     *
     * After all nodes and elements are parsed, buildConnectivity() is called
     * automatically on the mesh to finalise topology.
     *
     * @param filepath Path to the .msh file.
     * @param mesh     Output mesh to populate.
     * @throws std::runtime_error on file-not-found, bad format, or unsupported version.
     */
    void read(const std::string& filepath, UnstructuredMesh& mesh);

private:
    /** @brief Validates the MeshFormat section; throws on non-v2 or binary format. */
    static void parseFormat(std::istream& in);

    /**
     * @brief Parses PhysicalNames section.
     * @return Map from Gmsh physical-group tag to {dimension, name}.
     */
    [[nodiscard]] static std::unordered_map<int, std::pair<int, std::string>>
    parsePhysicalNames(std::istream& in);

    /**
     * @brief Parses the Nodes section, adding nodes to the mesh.
     * @param in            Input stream positioned after "$Nodes".
     * @param mesh          Mesh to populate with nodes.
     * @param gmshIdToNode  Output map from Gmsh 1-based node ID to NodeId.
     */
    static void parseNodes(std::istream& in,
                           UnstructuredMesh& mesh,
                           std::unordered_map<int, NodeId>& gmshIdToNode);

    /**
     * @brief Parses the Elements section, creating cells and boundary patches.
     *
     * Line elements (type 1) populate boundaryEdgeToPatch.
     * Triangle (type 2) and quad (type 3) elements are added as cells.
     *
     * @param in                  Input stream positioned after "$Elements".
     * @param mesh                Mesh to populate.
     * @param physGroups          Physical group map from parsePhysicalNames().
     * @param gmshIdToNode        Node ID translation map from parseNodes().
     * @param boundaryEdgeToPatch Output map from canonical edge key to patch index.
     */
    static void parseElements(
        std::istream&                                               in,
        UnstructuredMesh&                                          mesh,
        const std::unordered_map<int, std::pair<int, std::string>>& physGroups,
        std::unordered_map<int, NodeId>&                           gmshIdToNode,
        std::unordered_map<std::pair<int, int>, int, PairHash>&   boundaryEdgeToPatch);
};
