#include "io/VTKWriter.hpp"
#include "solver/NavierStokesSolver.hpp"
#include "utils/Config.hpp"
#include "utils/Logger.hpp"

#include <filesystem>
#include <stdexcept>
#include <string>

/**
 * @brief Entry point for the LG-Flow solver.
 *
 * Usage: lgflow <config-file>
 *
 * Loads the configuration, constructs the solver, runs the SIMPLE loop,
 * and writes VTK output on completion.
 */
int main(int argc, char* argv[])
{
    Logger& log = Logger::get();
    log.info("LG-Flow - 2D incompressible Navier-Stokes solver (stub)");

    if (argc < 2)
    {
        log.error("Usage: lgflow <config-file>");
        return 1;
    }

    const std::string configPath = argv[1];

    // Validate that the config file path exists before attempting to parse it.
    if (!std::filesystem::exists(configPath))
    {
        log.error("Config file not found: '" + configPath + "'");
        return 1;
    }

    try
    {
        // Load configuration
        Config config;
        config.load(configPath);
        log.info("Config loaded from: " + configPath);

        // Construct and initialise solver
        NavierStokesSolver solver(config);
        solver.initialize();

        // Run SIMPLE loop
        const int maxIter = config.get<int>("solver.max_iter", 1000);
        solver.run(maxIter);

        log.info("Solver finished. Residual = "
                 + std::to_string(solver.residual()));

        // Prepare output directory: create it if it does not exist.
        const std::string outDir = config.get<std::string>("output.dir", "output");
        std::error_code ec;
        std::filesystem::create_directories(outDir, ec);
        if (ec)
        {
            log.error("Cannot create output directory '" + outDir
                      + "': " + ec.message());
            return 1;
        }

        // Write VTK output
        VTKWriter writer;
        writer.write(outDir + "/result.vtu",
                     solver.velocity().mesh(),
                     solver.pressure(),
                     solver.velocity());

        log.info("VTK output written to " + outDir + "/result.vtu");
    }
    catch (const std::exception& e)
    {
        Logger::get().error(std::string("Fatal error: ") + e.what());
        return 1;
    }

    return 0;
}
