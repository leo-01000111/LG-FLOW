#include "io/VTKWriter.hpp"
#include "core/Field.hpp"
#include "core/Mesh.hpp"

#include <Eigen/Dense>
#include <filesystem>
#include <fstream>
#include <gtest/gtest.h>
#include <sstream>
#include <stdexcept>
#include <string>

// ── Helpers ───────────────────────────────────────────────────────────────────

namespace
{

/// Reads the entire contents of a file into a string.
std::string readFile(const std::string& path)
{
    std::ifstream in(path);
    std::ostringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

/// Returns a unique temp-file path for a given test name.
std::string tmpVtu(const std::string& name)
{
    namespace fs = std::filesystem;
    return (fs::temp_directory_path() / ("lgflow_vtktest_" + name + ".vtu"))
        .string();
}

}  // namespace

// ── File creation ──────────────────────────────────────────────────────────────

TEST(VTKWriter, Write_CreatesFile)
{
    Mesh mesh;
    mesh.load(4, 4, 1.0, 1.0);
    Field<double>          p(mesh, 0.0);
    Field<Eigen::Vector2d> v(mesh, Eigen::Vector2d::Zero());

    const std::string path = tmpVtu("creates_file");

    VTKWriter writer;
    writer.write(path, mesh, p, v);

    EXPECT_TRUE(std::filesystem::exists(path));
    std::filesystem::remove(path);
}

// ── XML header present ────────────────────────────────────────────────────────

TEST(VTKWriter, Write_ContainsXmlDeclaration)
{
    Mesh mesh;
    mesh.load(3, 3, 1.0, 1.0);
    Field<double>          p(mesh, 1.5);
    Field<Eigen::Vector2d> v(mesh, Eigen::Vector2d(1.0, 2.0));

    const std::string path = tmpVtu("xml_decl");
    VTKWriter().write(path, mesh, p, v);

    const std::string content = readFile(path);
    std::filesystem::remove(path);

    EXPECT_NE(content.find("<?xml"), std::string::npos);
    EXPECT_NE(content.find("VTKFile"), std::string::npos);
}

// ── UnstructuredGrid element ──────────────────────────────────────────────────

TEST(VTKWriter, Write_ContainsUnstructuredGrid)
{
    Mesh mesh;
    mesh.load(2, 2, 1.0, 1.0);
    Field<double>          p(mesh, 0.0);
    Field<Eigen::Vector2d> v(mesh, Eigen::Vector2d::Zero());

    const std::string path = tmpVtu("unstruc");
    VTKWriter().write(path, mesh, p, v);

    const std::string content = readFile(path);
    std::filesystem::remove(path);

    EXPECT_NE(content.find("UnstructuredGrid"), std::string::npos);
}

// ── Correct NumberOfPoints and NumberOfCells ──────────────────────────────────

TEST(VTKWriter, Write_NumberOfPointsAndCellsCorrect)
{
    const int Nx = 5, Ny = 3;
    const int Nc = Nx * Ny;  // 15

    Mesh mesh;
    mesh.load(Nx, Ny, 1.0, 1.0);
    Field<double>          p(mesh, 0.0);
    Field<Eigen::Vector2d> v(mesh, Eigen::Vector2d::Zero());

    const std::string path = tmpVtu("npts_ncells");
    VTKWriter().write(path, mesh, p, v);

    const std::string content = readFile(path);
    std::filesystem::remove(path);

    const std::string expectPts   = "NumberOfPoints=\"" + std::to_string(Nc) + "\"";
    const std::string expectCells = "NumberOfCells=\""  + std::to_string(Nc) + "\"";
    EXPECT_NE(content.find(expectPts),   std::string::npos) << "Missing: " << expectPts;
    EXPECT_NE(content.find(expectCells), std::string::npos) << "Missing: " << expectCells;
}

// ── Pressure array name ───────────────────────────────────────────────────────

TEST(VTKWriter, Write_PressureArrayNamePresent)
{
    Mesh mesh;
    mesh.load(2, 2, 1.0, 1.0);
    Field<double>          p(mesh, 3.14);
    Field<Eigen::Vector2d> v(mesh, Eigen::Vector2d::Zero());

    const std::string path = tmpVtu("pres_name");
    VTKWriter().write(path, mesh, p, v);

    const std::string content = readFile(path);
    std::filesystem::remove(path);

    EXPECT_NE(content.find("Name=\"pressure\""), std::string::npos);
}

// ── Velocity array name ───────────────────────────────────────────────────────

TEST(VTKWriter, Write_VelocityArrayNamePresent)
{
    Mesh mesh;
    mesh.load(2, 2, 1.0, 1.0);
    Field<double>          p(mesh, 0.0);
    Field<Eigen::Vector2d> v(mesh, Eigen::Vector2d(1.0, 0.5));

    const std::string path = tmpVtu("vel_name");
    VTKWriter().write(path, mesh, p, v);

    const std::string content = readFile(path);
    std::filesystem::remove(path);

    EXPECT_NE(content.find("Name=\"velocity\""), std::string::npos);
}

// ── Size mismatch validation ──────────────────────────────────────────────────

TEST(VTKWriter, Write_PressureSizeMismatch_Throws)
{
    Mesh mesh1, mesh2;
    mesh1.load(4, 4, 1.0, 1.0);  // 16 cells
    mesh2.load(6, 6, 1.0, 1.0);  // 36 cells

    Field<double>          pBig(mesh2, 0.0);             // size 36 — mismatch
    Field<Eigen::Vector2d> v(mesh1, Eigen::Vector2d::Zero());

    const std::string path = tmpVtu("pres_mismatch");
    VTKWriter writer;
    EXPECT_THROW(writer.write(path, mesh1, pBig, v), std::invalid_argument);
    std::filesystem::remove(path);  // cleanup if file was somehow created
}

TEST(VTKWriter, Write_VelocitySizeMismatch_Throws)
{
    Mesh mesh1, mesh2;
    mesh1.load(4, 4, 1.0, 1.0);  // 16 cells
    mesh2.load(6, 6, 1.0, 1.0);  // 36 cells

    Field<double>          p(mesh1, 0.0);
    Field<Eigen::Vector2d> vBig(mesh2, Eigen::Vector2d::Zero());  // size 36 — mismatch

    const std::string path = tmpVtu("vel_mismatch");
    VTKWriter writer;
    EXPECT_THROW(writer.write(path, mesh1, p, vBig), std::invalid_argument);
    std::filesystem::remove(path);
}

// ── Non-existent output directory ────────────────────────────────────────────

TEST(VTKWriter, Write_NonExistentDirectory_Throws)
{
    Mesh mesh;
    mesh.load(2, 2, 1.0, 1.0);
    Field<double>          p(mesh, 0.0);
    Field<Eigen::Vector2d> v(mesh, Eigen::Vector2d::Zero());

    // Path whose parent directory does not exist — ofstream will fail to open.
    namespace fs = std::filesystem;
    const std::string badPath =
        (fs::temp_directory_path() / "lgflow_no_such_subdir_xyz" / "test.vtu")
        .string();

    VTKWriter writer;
    EXPECT_THROW(writer.write(badPath, mesh, p, v), std::runtime_error);
}
