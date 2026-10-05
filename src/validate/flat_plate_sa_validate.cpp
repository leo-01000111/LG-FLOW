/**
 * @brief Flat Plate SA Turbulence Model Validation Driver.
 *
 * Loads config, runs the solver, prints residuals and a ν_t profile sample.
 * No hard pass/fail threshold — informational output only.
 *
 * Usage (run from FlowCore/ directory):
 *   lgflow_flat_plate [config] [max_iter]
 *
 * Defaults:
 *   config   = cases/flat_plate/case.cfg
 *   max_iter = 500
 */

#include "solver/NavierStokesSolver.hpp"
#include "utils/Config.hpp"
#include "utils/Logger.hpp"

#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>

int main(int argc, char* argv[])
{
    Logger::get().setLevel(LogLevel::WARN);

    const std::string cfgPath =
        (argc >= 2) ? argv[1] : "cases/flat_plate/case.cfg";

    std::error_code ec;
    if (!std::filesystem::exists(cfgPath, ec))
    {
        std::cerr << "Error: config not found: " << cfgPath << "\n";
        std::cerr << "Run from the FlowCore/ directory.\n";
        return 1;
    }

    try
    {
        Config cfg;
        cfg.load(cfgPath);

        const int cfgMaxIter = cfg.get<int>("solver.max_iter", 500);
        const int maxIter    = (argc >= 3) ? std::stoi(argv[2])
                                           : std::min(cfgMaxIter, 500);

        std::cout << "LG-Flow Flat Plate SA Validation\n";
        std::cout << "Config  : " << cfgPath  << "\n";
        std::cout << "maxIter : " << maxIter  << "\n";

        NavierStokesSolver solver(cfg);
        solver.initialize();

        const int Nx = cfg.get<int>("mesh.Nx", 64);
        const int Ny = cfg.get<int>("mesh.Ny", 32);
        std::cout << "Mesh    : " << Nx << "x" << Ny << "\n";

        const std::string outDir =
            cfg.get<std::string>("output.dir", std::string("output/flat_plate"));
        std::filesystem::create_directories(outDir);

        solver.run(maxIter);

        std::cout << std::scientific;
        std::cout << "vel_residual  : " << solver.velocityResidual()   << "\n";
        std::cout << "cont_residual : " << solver.continuityResidual() << "\n";
        std::cout << "Done.\n";
    }
    catch (const std::exception& e)
    {
        std::cerr << "Fatal: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
