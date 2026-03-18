#include "utils/Config.hpp"

#include <filesystem>
#include <fstream>
#include <gtest/gtest.h>
#include <stdexcept>
#include <string>

// ── Helper: write a string to a temporary file ────────────────────────────────

namespace
{
/// Writes content to a unique temp file and returns its path.
/// The caller is responsible for deleting the file when done.
std::string writeTempConfig(const std::string& content)
{
    auto tmp = std::filesystem::temp_directory_path()
               / ("lgflow_test_" + std::to_string(std::hash<std::string>{}(content)) + ".cfg");
    std::ofstream f(tmp);
    f << content;
    return tmp.string();
}
}  // namespace

// ── Load: basic parsing ───────────────────────────────────────────────────────

TEST(Config, Load_ValidKeyValue_ParsedCorrectly)
{
    std::string path = writeTempConfig("solver.dt = 0.05\nsolver.rho = 1.2\n");
    Config cfg;
    cfg.load(path);

    EXPECT_NEAR(cfg.get<double>("solver.dt"),  0.05, 1e-12);
    EXPECT_NEAR(cfg.get<double>("solver.rho"), 1.2,  1e-12);

    std::filesystem::remove(path);
}

TEST(Config, Load_CommentLines_Ignored)
{
    std::string path = writeTempConfig("# this is a comment\nsolver.dt = 0.1\n");
    Config cfg;
    cfg.load(path);

    EXPECT_NEAR(cfg.get<double>("solver.dt"), 0.1, 1e-12);

    std::filesystem::remove(path);
}

TEST(Config, Load_BlankLines_Ignored)
{
    std::string path = writeTempConfig("\n\nsolver.dt = 0.2\n\n");
    Config cfg;
    cfg.load(path);

    EXPECT_NEAR(cfg.get<double>("solver.dt"), 0.2, 1e-12);

    std::filesystem::remove(path);
}

TEST(Config, Load_LeadingTrailingWhitespace_Trimmed)
{
    std::string path = writeTempConfig("  solver.dt  =  0.3  \n");
    Config cfg;
    cfg.load(path);

    EXPECT_NEAR(cfg.get<double>("solver.dt"), 0.3, 1e-12);

    std::filesystem::remove(path);
}

TEST(Config, Load_IntValue_ParsedCorrectly)
{
    std::string path = writeTempConfig("mesh.Nx = 64\n");
    Config cfg;
    cfg.load(path);

    EXPECT_EQ(cfg.get<int>("mesh.Nx"), 64);

    std::filesystem::remove(path);
}

// ── Duplicate-key policy: last-wins ──────────────────────────────────────────

TEST(Config, Load_DuplicateKey_LastWins)
{
    // Policy: if a key appears more than once, the last value is used.
    std::string path = writeTempConfig("solver.dt = 0.01\nsolver.dt = 0.99\n");
    Config cfg;
    cfg.load(path);

    EXPECT_NEAR(cfg.get<double>("solver.dt"), 0.99, 1e-12);

    std::filesystem::remove(path);
}

// ── has() ─────────────────────────────────────────────────────────────────────

TEST(Config, Has_PresentKey_ReturnsTrue)
{
    std::string path = writeTempConfig("solver.dt = 0.01\n");
    Config cfg;
    cfg.load(path);

    EXPECT_TRUE(cfg.has("solver.dt"));

    std::filesystem::remove(path);
}

TEST(Config, Has_MissingKey_ReturnsFalse)
{
    Config cfg;  // empty, no file loaded
    EXPECT_FALSE(cfg.has("nonexistent.key"));
}

// ── get<T> with default ───────────────────────────────────────────────────────

TEST(Config, Get_MissingKeyWithDefault_ReturnsDefault)
{
    Config cfg;
    EXPECT_NEAR(cfg.get<double>("missing", 42.0), 42.0, 1e-12);
    EXPECT_EQ  (cfg.get<int>   ("missing", 7),    7);
}

TEST(Config, Get_PresentKeyIgnoresDefault)
{
    std::string path = writeTempConfig("mesh.Nx = 32\n");
    Config cfg;
    cfg.load(path);

    EXPECT_EQ(cfg.get<int>("mesh.Nx", 999), 32);

    std::filesystem::remove(path);
}

// ── Error cases ───────────────────────────────────────────────────────────────

TEST(Config, Load_MissingFile_Throws)
{
    Config cfg;
    EXPECT_THROW(cfg.load("/definitely/does/not/exist.cfg"), std::runtime_error);
}

TEST(Config, Load_MalformedLine_NoEquals_Throws)
{
    std::string path = writeTempConfig("this line has no equals sign\n");
    Config cfg;
    EXPECT_THROW(cfg.load(path), std::runtime_error);
    std::filesystem::remove(path);
}

TEST(Config, Load_EmptyKey_Throws)
{
    std::string path = writeTempConfig("  = value\n");
    Config cfg;
    EXPECT_THROW(cfg.load(path), std::runtime_error);
    std::filesystem::remove(path);
}

TEST(Config, Get_MissingKeyNoDefault_Throws)
{
    Config cfg;
    EXPECT_THROW(static_cast<void>(cfg.get<double>("nonexistent.key")), std::runtime_error);
}
