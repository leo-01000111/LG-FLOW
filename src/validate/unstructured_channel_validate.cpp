/**
 * @brief Unstructured Channel Flow Validation Driver.
 *
 * Loads the channel_unstruct case, runs up to 500 iterations, and prints
 * velocity and continuity residuals at the end.
 *
 * Usage (run from FlowCore/):
 *   lgflow_unstruct [config] [max_iter]
 *
 * Defaults:
 *   config   = cases/channel_unstruct/case.cfg
 *   max_iter = 500
 */

#include "solver/UnstructuredNavierStokesSolver.hpp"
#include "utils/Config.hpp"
#include "utils/Logger.hpp"

#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>

int main(int argc, char* argv[])
{
    Logger::get().setLevel(LogLevel::WARN);

    const std::string cfgPath =
        (argc >= 2) ? argv[1] : "cases/channel_unstruct/case.cfg";

    std::error_code ecExists;
    if (!std::filesystem::exists(cfgPath, ecExists))
    {
        std::cerr << "Error: config not found: " << cfgPath << "\n";
        std::cerr << "Run from the FlowCore/ directory.\n";
        return 1;
    }

    const int maxIter = (argc >= 3) ? std::stoi(argv[2]) : 500;

    try
    {
        Config cfg;
        cfg.load(cfgPath);

        std::cout << "LG-Flow Unstructured: Channel Flow Validation\n";
        std::cout << "Config  : " << cfgPath << "\n";
        std::cout << "maxIter : " << maxIter << "\n";

        UnstructuredNavierStokesSolver solver(cfg);
        solver.initialize();

        std::cout << "Mesh    : "
                  << solver.velocity().mesh().numCells() << " cells, "
                  << solver.velocity().mesh().numFaces() << " faces\n";

        solver.run(maxIter);

        std::cout << std::scientific;
        std::cout << "Results : vel_residual=" << solver.velocityResidual()
                  << "  cont_residual=" << solver.continuityResidual() << "\n";
    }
    catch (const std::exception& e)
    {
        std::cerr << "Fatal: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
