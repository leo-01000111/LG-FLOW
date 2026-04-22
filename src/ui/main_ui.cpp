#include "ResidualPlot.hpp"
#include "SidePanel.hpp"
#include "StatusBar.hpp"
#include "UIColors.hpp"
#include "UILayout.hpp"

#include "solver/NavierStokesSolver.hpp"
#include "utils/Config.hpp"
#include "utils/Logger.hpp"

#include <SFML/Graphics.hpp>
#include <SFML/Window/VideoMode.hpp>

#include <atomic>
#include <filesystem>
#include <mutex>
#include <optional>
#include <string>
#include <thread>
#include <vector>

// ── SolverState ───────────────────────────────────────────────────────────────
// Shared between the UI thread and the solver thread. All fields accessed
// from both threads are protected by mutex or are atomic.
struct SolverState
{
    std::mutex              mtx;
    std::vector<double>     velHistory;
    std::vector<double>     contHistory;
    double                  lastVel   = 0.0;
    double                  lastCont  = 0.0;
    int                     iteration = 0;
    bool                    converged = false;

    std::atomic<bool>       stopRequest{false};
    std::atomic<bool>       running{false};
};

// ── solverThread ──────────────────────────────────────────────────────────────
// Runs the SIMPLE loop in a background thread. Posts residuals to SolverState
// after every step so the UI can display live convergence.
static void solverThread(const std::string& configPath, SolverState* state)
{
    state->running.store(true);
    state->stopRequest.store(false);

    try {
        Config cfg;
        cfg.load(configPath);

        const int maxIter = cfg.get<int>("solver.max_iter", 3000);
        const double tol  = cfg.get<double>("solver.tolerance", 1e-5);

        NavierStokesSolver solver(cfg);
        solver.initialize();

        for (int iter = 0; iter < maxIter; ++iter) {
            if (state->stopRequest.load()) break;

            solver.step(cfg.get<double>("solver.dt", 0.01));

            const double vel  = solver.velocityResidual();
            const double cont = solver.continuityResidual();

            {
                std::lock_guard<std::mutex> lk(state->mtx);
                state->velHistory.push_back(vel);
                state->contHistory.push_back(cont);
                state->lastVel   = vel;
                state->lastCont  = cont;
                state->iteration = iter + 1;
            }

            if (vel < tol) {
                std::lock_guard<std::mutex> lk(state->mtx);
                state->converged = true;
                break;
            }
        }
    } catch (const std::exception& ex) {
        Logger::get().error(std::string("Solver error: ") + ex.what());
    }

    state->running.store(false);
}

// ── main ──────────────────────────────────────────────────────────────────────
int main()
{
    // Discover available cases
    const std::vector<std::pair<std::string,std::string>> cases = {
        {"Lid-Driven Cavity Re=100", "cases/lid_driven_cavity/case_validate.cfg"},
        {"Channel Flow Re=100",      "cases/channel_flow/case.cfg"},
    };

    std::vector<std::string> caseNames;
    for (const auto& [name, _] : cases)
        caseNames.push_back(name);

    // Window
    sf::ContextSettings settings;
    settings.antiAliasingLevel = 8;
    sf::RenderWindow window(
        sf::VideoMode({static_cast<unsigned>(UILayout::WindowW),
                       static_cast<unsigned>(UILayout::WindowH)}),
        "LG-Flow CFD",
        sf::Style::Default,
        sf::State::Windowed,
        settings);
    window.setFramerateLimit(60);

    // Font — attempt to load system font, fall back gracefully
    sf::Font font;
    const std::vector<std::string> fontPaths = {
        "C:/Windows/Fonts/consola.ttf",
        "C:/Windows/Fonts/arial.ttf",
        "C:/Windows/Fonts/segoeui.ttf",
    };
    bool fontLoaded = false;
    for (const auto& fp : fontPaths) {
        if (std::filesystem::exists(fp)) {
            if (font.openFromFile(fp)) { fontLoaded = true; break; }
        }
    }
    if (!fontLoaded) {
        // Use SFML default font if none found
        Logger::get().warn("No system font found; text may not render.");
    }

    // UI widgets
    SidePanel    sidePanel(font);
    StatusBar    statusBar(font);
    const sf::FloatRect plotBounds{
        {UILayout::PlotX, UILayout::PlotY},
        {UILayout::PlotW, UILayout::PlotH}};
    ResidualPlot plot(plotBounds, font);

    sidePanel.setCases(caseNames);

    // State
    SolverState  state;
    std::thread  solverThr;
    int          selectedCase = 0;
    bool         prevRunning  = false;

    // Wire callbacks
    sidePanel.onCaseSelected = [&](int idx) {
        selectedCase = idx;
        // Reset plot when case changes
        std::lock_guard<std::mutex> lk(state.mtx);
        state.velHistory.clear();
        state.contHistory.clear();
        state.lastVel   = 0.0;
        state.lastCont  = 0.0;
        state.iteration = 0;
        state.converged = false;
        plot.clear();
    };

    sidePanel.onRunPressed = [&]() {
        if (state.running.load()) return;
        {
            std::lock_guard<std::mutex> lk(state.mtx);
            state.velHistory.clear();
            state.contHistory.clear();
            state.lastVel   = 0.0;
            state.lastCont  = 0.0;
            state.iteration = 0;
            state.converged = false;
        }
        plot.clear();
        const std::string cfgPath = cases[static_cast<std::size_t>(selectedCase)].second;
        if (solverThr.joinable()) solverThr.join();
        solverThr = std::thread(solverThread, cfgPath, &state);
    };

    sidePanel.onStopPressed = [&]() {
        state.stopRequest.store(true);
    };

    sf::Clock clock;

    while (window.isOpen()) {
        const float dt = clock.restart().asSeconds();

        // Events
        while (const std::optional<sf::Event> event = window.pollEvent()) {
            if (event->is<sf::Event::Closed>()) {
                state.stopRequest.store(true);
                if (solverThr.joinable()) solverThr.join();
                window.close();
            }
            sidePanel.handleEvent(*event);
        }

        const bool running = state.running.load();
        sidePanel.setRunning(running);

        // Update plot with any new data
        {
            std::lock_guard<std::mutex> lk(state.mtx);
            const std::size_t plotted = [&] {
                // Count how many we've already pushed (approximate via deque size)
                // Simple approach: push everything since last frame.
                // We track last seen index.
                static std::size_t lastIdx = 0;
                const std::size_t sz = state.velHistory.size();
                for (std::size_t i = lastIdx; i < sz; ++i)
                    plot.push(state.velHistory[i], state.contHistory[i]);
                lastIdx = sz;
                return sz;
            }();
            (void)plotted;
        }

        // Read status for display
        double vel = 0, cont = 0;
        int iter = 0;
        bool converged = false;
        {
            std::lock_guard<std::mutex> lk(state.mtx);
            vel      = state.lastVel;
            cont     = state.lastCont;
            iter     = state.iteration;
            converged = state.converged;
        }
        statusBar.update(iter, vel, cont, running, converged);

        // Mouse pos for hover animation
        const sf::Vector2i mp = sf::Mouse::getPosition(window);
        const sf::Vector2f mousePos{static_cast<float>(mp.x),
                                    static_cast<float>(mp.y)};
        sidePanel.update(dt, mousePos);

        // Draw
        window.clear(UIColors::Background);
        sidePanel.draw(window);
        statusBar.draw(window);
        plot.draw(window);
        window.display();

        prevRunning = running;
    }

    state.stopRequest.store(true);
    if (solverThr.joinable()) solverThr.join();
    return 0;
}
