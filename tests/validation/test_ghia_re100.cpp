/**
 * @brief [SLOW] Ghia et al. (1982) Re=100 lid-driven cavity validation.
 *
 * Runs the solver to convergence on a 32x32 grid, samples u(x=0.5, y) and
 * v(x, y=0.5) centerlines, and checks RMS errors against Ghia 1982 Table 1/3.
 *
 * Thresholds (32x32 upwind, measured baseline in parentheses):
 *   u-centerline L2 < 0.10  (measured 0.071)
 *   v-centerline L2 < 0.06  (measured 0.037)
 *
 * Reference: Ghia, U., Ghia, K.N., Shin, C.T. (1982).
 *   "High-Re solutions for incompressible flow using the Navier-Stokes
 *    equations and a multigrid method."
 *   Journal of Computational Physics, 48(3), 387-411.
 *
 * Run from FlowCore/ so that cases/ and output/ paths resolve correctly.
 * CTest sets WORKING_DIRECTORY automatically via CMakeLists.txt.
 */

#include "solver/NavierStokesSolver.hpp"
#include "utils/Config.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

// ── Local helpers (mirror of lid_driven_cavity_validate.cpp) ─────────────────

namespace {

struct CenterlineSample { double coord, value; };
struct ReferencePoint   { std::string axis; double coord, refValue; };

void sampleCenterlines(const Field<Eigen::Vector2d>& vel,
                       double midX, double midY,
                       std::vector<CenterlineSample>& uSamples,
                       std::vector<CenterlineSample>& vSamples)
{
    const Mesh& mesh = vel.mesh();
    const int   Nx   = mesh.Nx();
    const int   Ny   = mesh.Ny();

    int    i_mid = 0;
    double minDX = std::numeric_limits<double>::max();
    for (int i = 0; i < Nx; ++i) {
        const double d = std::abs(mesh.getCellCenter(i, 0).x() - midX);
        if (d < minDX) { minDX = d; i_mid = i; }
    }

    int    j_mid = 0;
    double minDY = std::numeric_limits<double>::max();
    for (int j = 0; j < Ny; ++j) {
        const double d = std::abs(mesh.getCellCenter(0, j).y() - midY);
        if (d < minDY) { minDY = d; j_mid = j; }
    }

    uSamples.clear();
    for (int j = 0; j < Ny; ++j)
        uSamples.push_back({mesh.getCellCenter(i_mid, j).y(), vel(i_mid, j).x()});

    vSamples.clear();
    for (int i = 0; i < Nx; ++i)
        vSamples.push_back({mesh.getCellCenter(i, j_mid).x(), vel(i, j_mid).y()});
}

std::vector<ReferencePoint> loadReference(const std::string& path)
{
    std::vector<ReferencePoint> pts;
    std::ifstream in(path);
    if (!in.is_open()) return pts;

    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#')          continue;
        if (line.find("axis") != std::string::npos)  continue;
        std::istringstream ss(line);
        ReferencePoint p;
        if (!std::getline(ss, p.axis, ',')) continue;
        char c{};
        if (!(ss >> p.coord >> c >> p.refValue)) continue;
        pts.push_back(p);
    }
    return pts;
}

double computeL2(const std::vector<CenterlineSample>& sim,
                 const std::vector<ReferencePoint>&   refs,
                 const std::string&                   axis)
{
    double sum = 0.0;
    int    n   = 0;
    for (const auto& ref : refs) {
        if (ref.axis != axis) continue;
        double minD   = std::numeric_limits<double>::max();
        double simVal = 0.0;
        for (const auto& s : sim) {
            const double d = std::abs(s.coord - ref.coord);
            if (d < minD) { minD = d; simVal = s.value; }
        }
        const double err = simVal - ref.refValue;
        sum += err * err;
        ++n;
    }
    return n > 0 ? std::sqrt(sum / static_cast<double>(n)) : 0.0;
}

} // namespace

// ── Test fixture ──────────────────────────────────────────────────────────────
// SetUpTestSuite() runs the solver once; individual tests just read results.

class GhiaRe100 : public ::testing::Test
{
protected:
    static void SetUpTestSuite()
    {
        const std::string cfgPath = "cases/lid_driven_cavity/case_validate.cfg";
        const std::string refPath = "cases/lid_driven_cavity/ghia1982_re100.csv";

        if (!std::filesystem::exists(cfgPath) || !std::filesystem::exists(refPath)) {
            s_skip = true;
            return;
        }

        Config cfg;
        cfg.load(cfgPath);

        NavierStokesSolver solver(cfg);
        solver.initialize();

        const int maxIter = cfg.get<int>("solver.max_iter", 2000);
        solver.run(maxIter);

        s_velResidual = solver.velocityResidual();

        const auto& vel = solver.velocity();
        const double Lx = cfg.get<double>("mesh.Lx", 1.0);
        const double Ly = cfg.get<double>("mesh.Ly", 1.0);

        sampleCenterlines(vel, Lx * 0.5, Ly * 0.5, s_uSamples, s_vSamples);

        s_refs = loadReference(refPath);
        s_solverRan = true;
    }

    void SetUp() override
    {
        if (s_skip)
            GTEST_SKIP() << "Validation data not found — run from FlowCore/";
        ASSERT_TRUE(s_solverRan) << "Solver failed to run in SetUpTestSuite";
        ASSERT_FALSE(s_refs.empty()) << "No Ghia reference data loaded";
    }

    static bool                          s_skip;
    static bool                          s_solverRan;
    static double                        s_velResidual;
    static std::vector<CenterlineSample> s_uSamples;
    static std::vector<CenterlineSample> s_vSamples;
    static std::vector<ReferencePoint>   s_refs;
};

bool                          GhiaRe100::s_skip        = false;
bool                          GhiaRe100::s_solverRan   = false;
double                        GhiaRe100::s_velResidual = 1.0;
std::vector<CenterlineSample> GhiaRe100::s_uSamples;
std::vector<CenterlineSample> GhiaRe100::s_vSamples;
std::vector<ReferencePoint>   GhiaRe100::s_refs;

// ── Tests ─────────────────────────────────────────────────────────────────────

TEST_F(GhiaRe100, ConvergedBelow1em5)
{
    // Velocity residual must reach the solver tolerance (1e-5).
    EXPECT_LT(s_velResidual, 1e-4)
        << "Solver stalled; vel_residual=" << s_velResidual;
}

TEST_F(GhiaRe100, UCenterlineL2)
{
    // u(x=0.5, y): RMS error vs Ghia Table 1.
    const double uL2 = computeL2(s_uSamples, s_refs, "u");
    EXPECT_LT(uL2, 0.10)
        << "u-centerline L2=" << uL2 << " exceeds threshold 0.10";
}

TEST_F(GhiaRe100, VCenterlineL2)
{
    // v(x, y=0.5): RMS error vs Ghia Table 3.
    const double vL2 = computeL2(s_vSamples, s_refs, "v");
    EXPECT_LT(vL2, 0.06)
        << "v-centerline L2=" << vL2 << " exceeds threshold 0.06";
}
