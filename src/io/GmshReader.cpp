#include "io/GmshReader.hpp"

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// Gmsh .msh format v2 reader.
// Reference: http://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format-version-2

void GmshReader::read(const std::string& filepath, UnstructuredMesh& mesh)
{
    std::ifstream in(filepath);
    if (!in.is_open())
        throw std::runtime_error("GmshReader::read: cannot open file '" + filepath + "'");

    std::unordered_map<int, std::pair<int, std::string>> physGroups;
    std::unordered_map<int, NodeId>                      gmshIdToNode;
    std::unordered_map<std::pair<int, int>, int, PairHash> boundaryEdgeToPatch;

    std::string line;
    while (std::getline(in, line))
    {
        if (line == "$MeshFormat")
            parseFormat(in);
        else if (line == "$PhysicalNames")
            physGroups = parsePhysicalNames(in);
        else if (line == "$Nodes")
            parseNodes(in, mesh, gmshIdToNode);
        else if (line == "$Elements")
            parseElements(in, mesh, physGroups, gmshIdToNode, boundaryEdgeToPatch);
    }

    if (mesh.numCells() == 0)
        throw std::runtime_error(
            "GmshReader::read: no 2D elements found in '" + filepath + "'");

    mesh.buildConnectivity(boundaryEdgeToPatch);
}

void GmshReader::parseFormat(std::istream& in)
{
    std::string line;
    if (!std::getline(in, line))
        throw std::runtime_error("GmshReader: unexpected EOF reading MeshFormat");

    std::istringstream ss(line);
    double version{};
    int fileType{};
    int dataSize{};
    if (!(ss >> version >> fileType >> dataSize))
        throw std::runtime_error("GmshReader: malformed MeshFormat line: '" + line + "'");

    if (version < 2.0 || version >= 3.0)
        throw std::runtime_error(
            "GmshReader: unsupported Gmsh format version "
            + std::to_string(version) + " (expected 2.x)");

    if (fileType != 0)
        throw std::runtime_error(
            "GmshReader: binary Gmsh files are not supported (file-type=" +
            std::to_string(fileType) + "); use ASCII format");
}

std::unordered_map<int, std::pair<int, std::string>>
GmshReader::parsePhysicalNames(std::istream& in)
{
    std::unordered_map<int, std::pair<int, std::string>> groups;

    std::string line;
    if (!std::getline(in, line))
        throw std::runtime_error("GmshReader: unexpected EOF reading PhysicalNames count");

    std::istringstream countSS(line);
    int count{};
    if (!(countSS >> count))
        throw std::runtime_error("GmshReader: bad PhysicalNames count line: '" + line + "'");

    for (int i = 0; i < count; ++i)
    {
        if (!std::getline(in, line))
            throw std::runtime_error("GmshReader: unexpected EOF in PhysicalNames");

        std::istringstream ss(line);
        int  dim{};
        int  tag{};
        std::string name;
        if (!(ss >> dim >> tag))
            throw std::runtime_error("GmshReader: bad PhysicalNames entry: '" + line + "'");

        // Name is quoted; read the rest of the line and strip quotes.
        if (std::getline(ss >> std::ws, name))
        {
            if (!name.empty() && name.front() == '"')
                name.erase(name.begin());
            if (!name.empty() && name.back() == '"')
                name.pop_back();
        }

        groups[tag] = {dim, name};
    }

    return groups;
}

void GmshReader::parseNodes(std::istream&                     in,
                             UnstructuredMesh&                 mesh,
                             std::unordered_map<int, NodeId>&  gmshIdToNode)
{
    std::string line;
    if (!std::getline(in, line))
        throw std::runtime_error("GmshReader: unexpected EOF reading Nodes count");

    std::istringstream countSS(line);
    int count{};
    if (!(countSS >> count))
        throw std::runtime_error("GmshReader: bad Nodes count line: '" + line + "'");

    gmshIdToNode.reserve(static_cast<std::size_t>(count));

    for (int i = 0; i < count; ++i)
    {
        if (!std::getline(in, line))
            throw std::runtime_error("GmshReader: unexpected EOF in Nodes");

        std::istringstream ss(line);
        int    gmshId{};
        double x{};
        double y{};
        double z{};

        if (!(ss >> gmshId >> x >> y >> z))
            throw std::runtime_error("GmshReader: bad node line: '" + line + "'");

        MeshNode n;
        n.position = Eigen::Vector2d(x, y);
        // z is discarded (2D solver)

        const NodeId nid = mesh.addNode(n);
        gmshIdToNode[gmshId] = nid;
    }
}

void GmshReader::parseElements(
    std::istream&                                                in,
    UnstructuredMesh&                                           mesh,
    const std::unordered_map<int, std::pair<int, std::string>>& physGroups,
    std::unordered_map<int, NodeId>&                            gmshIdToNode,
    std::unordered_map<std::pair<int, int>, int, PairHash>&    boundaryEdgeToPatch)
{
    std::string line;
    if (!std::getline(in, line))
        throw std::runtime_error("GmshReader: unexpected EOF reading Elements count");

    std::istringstream countSS(line);
    int count{};
    if (!(countSS >> count))
        throw std::runtime_error("GmshReader: bad Elements count line: '" + line + "'");

    for (int i = 0; i < count; ++i)
    {
        if (!std::getline(in, line))
            throw std::runtime_error("GmshReader: unexpected EOF in Elements");

        std::istringstream ss(line);
        int elemId{};
        int elemType{};
        int nTags{};

        if (!(ss >> elemId >> elemType >> nTags))
            throw std::runtime_error("GmshReader: bad element header: '" + line + "'");

        // Read tags; first tag is always the physical group ID.
        int physGroupTag = 0;
        for (int t = 0; t < nTags; ++t)
        {
            int tag{};
            if (!(ss >> tag))
                throw std::runtime_error("GmshReader: expected tag in element line: '" + line + "'");
            if (t == 0)
                physGroupTag = tag;
        }

        // Determine number of nodes for this element type.
        int numElemNodes = 0;
        switch (elemType)
        {
            case 1:  numElemNodes = 2; break;  // line
            case 2:  numElemNodes = 3; break;  // triangle
            case 3:  numElemNodes = 4; break;  // quad
            case 15: numElemNodes = 1; break;  // point (ignored)
            default:
                // Skip unknown element types silently by reading their node IDs.
                // We don't know numElemNodes, so we can't skip cleanly —
                // treat as a parse error for robustness.
                throw std::runtime_error(
                    "GmshReader: unsupported element type " + std::to_string(elemType)
                    + " (only types 1=line, 2=triangle, 3=quad supported)");
        }

        std::vector<int> gmshNodeIds(static_cast<std::size_t>(numElemNodes));
        for (int k = 0; k < numElemNodes; ++k)
        {
            if (!(ss >> gmshNodeIds[static_cast<std::size_t>(k)]))
                throw std::runtime_error(
                    "GmshReader: missing node ID in element line: '" + line + "'");
        }

        if (elemType == 15)
            continue;  // point — ignore

        if (elemType == 1)
        {
            // Boundary edge — register patch and add to boundary map.
            const auto pgIt = physGroups.find(physGroupTag);
            if (pgIt == physGroups.end())
                continue;  // no physical group — skip

            const auto& [dim, name] = pgIt->second;
            if (dim != 1)
                continue;  // not a boundary group

            const int patchIdx = mesh.addPatch(name);

            const int n0 = toInt(gmshIdToNode.at(gmshNodeIds[0]));
            const int n1 = toInt(gmshIdToNode.at(gmshNodeIds[1]));
            const std::pair<int, int> key = (n0 < n1)
                ? std::pair<int, int>{n0, n1}
                : std::pair<int, int>{n1, n0};

            boundaryEdgeToPatch[key] = patchIdx;
        }
        else
        {
            // 2D element (triangle or quad) — add as cell.
            std::vector<NodeId> cellNodes(static_cast<std::size_t>(numElemNodes));
            for (int k = 0; k < numElemNodes; ++k)
                cellNodes[static_cast<std::size_t>(k)] = gmshIdToNode.at(gmshNodeIds[static_cast<std::size_t>(k)]);

            mesh.addCell(std::span<const NodeId>{cellNodes.data(), cellNodes.size()});
        }
    }
}
