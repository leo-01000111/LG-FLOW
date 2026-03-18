/**
 * @brief Lid-Driven Cavity Validation Executable (Phase 7).
 *
 * Runs the lid-driven cavity case, samples centerline velocities, writes a
 * CSV of sampled data, and compares against Ghia et al. (1982) reference if
 * the reference file is present.
 *
 * Usage (run from FlowCore/ directory):
 *   lgflow_validate_lid [config] [max_iter]
 *
 * Defaults:
 *   config   = cases/lid_driven_cavity/case_validate.cfg  (32x32, dt=0.001)
 *   max_iter = min(solver.max_iter from config, 500)
 *
 * Outputs:
 *   output/lid_driven_cavity/centerline.csv
 *   output/lid_driven_cavity/history.csv
 *   output/lid_driven_cavity/iter_*.vtu  (VTK snapshots)
 *
 * Sampling policy:
 *   u-centerline: u-velocity at the column i whose centre x is closest to 0.5.
 *                 Reported at each cell centre y in [0, Ly].
 *   v-centerline: v-velocity at the row j whose centre y is closest to 0.5.
 *                 Reported at each cell centre x in [0, Lx].
 *
 * Benchmark comparison:
 *   Reference file: cases/lid_driven_cavity/ghia1982_re100.csv
 *   Matching: nearest-neighbour by coord value.
 *   Metrics: L2 (sqrt(mean squared error)) and Linf (max absolute error).
 *   No hard threshold is enforced — metrics are reported informally.
 */

#include "solver/NavierStokesSolver.hpp"
#include "utils/Config.hpp"
#include "utils/Logger.hpp"

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// ── Data types ────────────────────────────────────────────────────────────────

struct CenterlineSample
{
    double coord;  ///< y (for u-line) or x (for v-line) [0, 1]
    double value;  ///< velocity component value [m/s]
};

struct ReferencePoint
{
    std::string axis;      ///< "u" or "v"
    double      coord;     ///< normalised coordinate
    double      refValue;  ///< reference velocity [m/s]
};

struct ErrorMetrics
{
    double l2{0.0};   ///< sqrt(mean squared error)
    double linf{0.0}; ///< maximum absolute error
    int    n{0};      ///< number of comparison points
};

// ── Centerline sampling ───────────────────────────────────────────────────────

/// Samples u(x≈midX, y) and v(x, y≈midY) from the post-run velocity field.
/// The column/row closest to midX or midY is selected (nearest-neighbour).
/// @param midX  x-coordinate of the vertical centerline (typically Lx/2)
/// @param midY  y-coordinate of the horizontal centerline (typically Ly/2)
void sampleCenterlines(const Field<Eigen::Vector2d>& vel,
                       double midX,
                       double midY,
                       std::vector<CenterlineSample>& uSamples,
                       std::vector<CenterlineSample>& vSamples)
{
    const Mesh& mesh = vel.mesh();
    const int   Nx   = mesh.Nx();
    const int   Ny   = mesh.Ny();

    // Find column i_mid: centre x closest to midX
    int    i_mid    = 0;
    double minDistX = std::numeric_limits<double>::max();
    for (int i = 0; i < Nx; ++i)
    {
        const double dist = std::abs(mesh.getCellCenter(i, 0).x() - midX);
        if (dist < minDistX) { minDistX = dist; i_mid = i; }
    }

    // Find row j_mid: centre y closest to midY
    int    j_mid    = 0;
    double minDistY = std::numeric_limits<double>::max();
    for (int j = 0; j < Ny; ++j)
    {
        const double dist = std::abs(mesh.getCellCenter(0, j).y() - midY);
        if (dist < minDistY) { minDistY = dist; j_mid = j; }
    }

    // u-centerline: u-velocity at column i_mid, for every row j
    uSamples.clear();
    uSamples.reserve(static_cast<std::size_t>(Ny));
    for (int j = 0; j < Ny; ++j)
        uSamples.push_back({mesh.getCellCenter(i_mid, j).y(), vel(i_mid, j).x()});

    // v-centerline: v-velocity at row j_mid, for every column i
    vSamples.clear();
    vSamples.reserve(static_cast<std::size_t>(Nx));
    for (int i = 0; i < Nx; ++i)
        vSamples.push_back({mesh.getCellCenter(i, j_mid).x(), vel(i, j_mid).y()});
}

// ── CSV I/O ───────────────────────────────────────────────────────────────────

/// Writes centerline samples to a CSV file.
/// Format: axis,coord,value  (one row per sample; u-samples before v-samples).
void writeCenterlineCSV(const std::string&                   path,
                        const std::vector<CenterlineSample>& uSamples,
                        const std::vector<CenterlineSample>& vSamples)
{
    std::ofstream out(path);
    if (!out.is_open())
        throw std::runtime_error("Cannot write centerline CSV: " + path);

    out << std::scientific << std::setprecision(8);
    out << "axis,coord,value\n";
    for (const auto& s : uSamples)
        out << "u," << s.coord << "," << s.value << "\n";
    for (const auto& s : vSamples)
        out << "v," << s.coord << "," << s.value << "\n";
}

/// Loads reference points from CSV.  Skips comment lines (starting with '#')
/// and the header line (containing "axis").  Returns an empty vector if the
/// file is not found.
std::vector<ReferencePoint> loadReference(const std::string& path)
{
    std::vector<ReferencePoint> pts;
    std::ifstream in(path);
    if (!in.is_open())
        return pts;

    std::string line;
    while (std::getline(in, line))
    {
        if (line.empty() || line[0] == '#')
            continue;
        if (line.find("axis") != std::string::npos)
            continue;  // header row

        std::istringstream ss(line);
        ReferencePoint     p;
        if (!std::getline(ss, p.axis, ','))
            continue;
        char comma{};
        if (!(ss >> p.coord >> comma >> p.refValue))
            continue;
        pts.push_back(p);
    }
    return pts;
}

// ── Error metrics ─────────────────────────────────────────────────────────────

/// Computes L2 and Linf errors between simulation samples and reference points
/// for a given axis ("u" or "v").  Nearest-neighbour matching by coord.
ErrorMetrics computeErrors(const std::vector<CenterlineSample>& sim,
                           const std::vector<ReferencePoint>&   refs,
                           const std::string&                   axis)
{
    ErrorMetrics m;
    for (const auto& ref : refs)
    {
        if (ref.axis != axis)
            continue;

        // Find nearest simulation sample by coordinate distance.
        double minDist = std::numeric_limits<double>::max();
        double simVal  = 0.0;
        for (const auto& s : sim)
        {
            const double d = std::abs(s.coord - ref.coord);
            if (d < minDist) { minDist = d; simVal = s.value; }
        }

        const double err = std::abs(simVal - ref.refValue);
        m.l2  += err * err;
        m.linf = std::max(m.linf, err);
        ++m.n;
    }
    if (m.n > 0)
        m.l2 = std::sqrt(m.l2 / static_cast<double>(m.n));  // RMS
    return m;
}

// ── Entry point ───────────────────────────────────────────────────────────────

int main(int argc, char* argv[])
{
    // Suppress INFO-level Logger output so only our structured stdout is visible.
    Logger::get().setLevel(LogLevel::WARN);

    const std::string cfgPath =
        (argc >= 2) ? argv[1] : "cases/lid_driven_cavity/case_validate.cfg";

    // Verify config file exists.
    std::error_code ecExists;
    if (!std::filesystem::exists(cfgPath, ecExists))
    {
        std::cerr << "Error: config not found: " << cfgPath << "\n";
        std::cerr << "Run from the FlowCore/ directory.\n";
        return 1;
    }

    try
    {
        Config cfg;
        cfg.load(cfgPath);

        // Determine iteration count: CLI arg overrides config, both capped
        // at 200 by default for CI speed. Pass a large CLI arg for full runs.
        const int cfgMaxIter = cfg.get<int>("solver.max_iter", 500);
        const int maxIter    = (argc >= 3)
                                ? std::stoi(argv[2])
                                : std::min(cfgMaxIter, 500);

        // Print header
        std::cout << "LG-Flow Validation: Lid-Driven Cavity\n";
        std::cout << "Config  : " << cfgPath << "\n";
        std::cout << "maxIter : " << maxIter << "\n";

        // Build and run solver
        NavierStokesSolver solver(cfg);
        solver.initialize();

        const int    Nx = cfg.get<int>   ("mesh.Nx", 16);
        const int    Ny = cfg.get<int>   ("mesh.Ny", 16);
        const double Lx = cfg.get<double>("mesh.Lx", 1.0);
        const double Ly = cfg.get<double>("mesh.Ly", 1.0);
        std::cout << "Mesh    : " << Nx << "x" << Ny << "\n";

        const std::string outDir =
            cfg.get<std::string>("output.dir", "output/lid_driven_cavity");
        std::filesystem::create_directories(outDir);

        solver.run(maxIter);

        // Count actual iterations from history.csv (written by run()).
        int actualIter = 0;
        {
            std::ifstream hist(outDir + "/history.csv");
            std::string   line;
            while (std::getline(hist, line))
                ++actualIter;
            if (actualIter > 0)
                --actualIter;  // subtract header row
        }

        std::cout << std::scientific << std::setprecision(3);
        std::cout << "Results : vel_residual=" << solver.velocityResidual()
                  << "  cont_residual=" << solver.continuityResidual()
                  << "  iterations=" << actualIter << "\n";

        // Sample centerlines from final velocity field.
        // Midpoints are at Lx/2 and Ly/2 to generalise beyond unit-square domains.
        const auto& vel = solver.velocity();
        std::vector<CenterlineSample> uSamples, vSamples;
        sampleCenterlines(vel, Lx * 0.5, Ly * 0.5, uSamples, vSamples);

        // Write centerline CSV
        const std::string csvPath = outDir + "/centerline.csv";
        writeCenterlineCSV(csvPath, uSamples, vSamples);
        std::cout << "CSV     : " << csvPath
                  << "  (" << uSamples.size() << " u-points, "
                  << vSamples.size() << " v-points)\n";

        // Load reference data and compute error metrics
        const std::string refPath =
            "cases/lid_driven_cavity/ghia1982_re100.csv";
        const auto refs = loadReference(refPath);

        ErrorMetrics uErr, vErr;

        if (refs.empty())
        {
            std::cout << "Ref     : not found (" << refPath << ")\n";
        }
        else
        {
            uErr = computeErrors(uSamples, refs, "u");
            vErr = computeErrors(vSamples, refs, "v");

            std::cout << "Ref     : " << refPath << "\n";
            std::cout << "u-line  : L2=" << uErr.l2
                      << "  Linf=" << uErr.linf
                      << "  (n=" << uErr.n << ")\n";
            std::cout << "v-line  : L2=" << vErr.l2
                      << "  Linf=" << vErr.linf
                      << "  (n=" << vErr.n << ")\n";
            std::cout << "[INFO] Metrics are informational — no pass/fail"
                         " threshold enforced.\n";
        }

        // Write validation_metrics.csv — machine-readable summary.
        {
            const std::string metricsPath = outDir + "/validation_metrics.csv";
            std::ofstream     mout(metricsPath);
            if (!mout.is_open())
                throw std::runtime_error(
                    "Cannot write validation_metrics.csv: " + metricsPath);

            mout << std::scientific << std::setprecision(8);
            mout << "metric,value\n";
            mout << "final_vel_residual,"  << solver.velocityResidual()    << "\n";
            mout << "final_cont_residual," << solver.continuityResidual()  << "\n";
            mout << "u_l2,"               << uErr.l2                      << "\n";
            mout << "u_linf,"             << uErr.linf                     << "\n";
            mout << "v_l2,"               << vErr.l2                      << "\n";
            mout << "v_linf,"             << vErr.linf                     << "\n";
            mout << "iterations,"         << actualIter                    << "\n";
            std::cout << "Metrics : " << metricsPath << "\n";
        }
    }
    catch (const std::exception& e)
    {
        std::cerr << "Fatal: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
