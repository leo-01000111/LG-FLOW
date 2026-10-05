#include "Colormap.hpp"
#include "FieldView.hpp"
#include "ResidualPlot.hpp"
#include "SidePanel.hpp"
#include "StatusBar.hpp"
#include "UIColors.hpp"
#include "UILayout.hpp"

#include "solver/NavierStokesSolver.hpp"
#include "utils/Config.hpp"
#include "utils/Logger.hpp"

#include <SFML/Graphics.hpp>
#include <SFML/Graphics/View.hpp>
#include <SFML/Window/VideoMode.hpp>

#include <atomic>
#include <filesystem>
#include <mutex>
#include <optional>
#include <string>
#include <thread>
#include <vector>

// ── CaseOverrides ─────────────────────────────────────────────────────────────
struct CaseOverrides {
    int    Nx      = 0;    // 0 = use config default
    int    Ny      = 0;
    double Re      = 0.0;  // 0.0 = use config default
    int    maxIter = 0;
};

// ── SolverState ───────────────────────────────────────────────────────────────
struct SolverState
{
    std::mutex              mtx;
    std::vector<double>     velHistory;
    std::vector<double>     contHistory;
    double                  lastVel   = 0.0;
    double                  lastCont  = 0.0;
    int                     iteration = 0;
    bool                    converged = false;

    FieldSnapshot           fieldSnap;
    bool                    snapReady = false;

    std::atomic<bool>       stopRequest{false};
    std::atomic<bool>       running{false};
};

// ── solverThread ──────────────────────────────────────────────────────────────
static void solverThread(const std::string& configPath, SolverState* state,
                         CaseOverrides ov)
{
    state->running.store(true);
    state->stopRequest.store(false);

    try {
        Config cfg;
        cfg.load(configPath);

        // Apply UI overrides before constructing the solver
        if (ov.Nx > 0)      cfg.set("mesh.Nx",        ov.Nx);
        if (ov.Ny > 0)      cfg.set("mesh.Ny",        ov.Ny);
        if (ov.Re > 0.0) {
            // Re = U_ref * Ly / nu; U_ref = 1.0 for all current cases
            const double Ly = cfg.get<double>("mesh.Ly", 1.0);
            cfg.set("solver.nu", Ly / ov.Re);
        }
        if (ov.maxIter > 0) cfg.set("solver.max_iter", ov.maxIter);

        const int    maxIter = cfg.get<int>("solver.max_iter", 3000);
        const double tol     = cfg.get<double>("solver.tolerance", 1e-5);
        const int    Nx      = cfg.get<int>("mesh.Nx", 16);
        const int    Ny      = cfg.get<int>("mesh.Ny", 16);
        const float  Lx      = static_cast<float>(cfg.get<double>("mesh.Lx", 1.0));
        const float  Ly      = static_cast<float>(cfg.get<double>("mesh.Ly", 1.0));

        NavierStokesSolver solver(cfg);
        solver.initialize();

        for (int iter = 0; iter < maxIter; ++iter) {
            if (state->stopRequest.load()) break;

            solver.step(cfg.get<double>("solver.dt", 0.01));

            const double vel  = solver.velocityResidual();
            const double cont = solver.continuityResidual();

            // Post residuals every step
            {
                std::lock_guard<std::mutex> lk(state->mtx);
                state->velHistory.push_back(vel);
                state->contHistory.push_back(cont);
                state->lastVel   = vel;
                state->lastCont  = cont;
                state->iteration = iter + 1;
            }

            // Post field snapshot every 25 iterations or on convergence
            const bool doSnap = (iter % 25 == 0) || (vel < tol);
            if (doSnap) {
                const auto& p = solver.pressure();
                const auto& u = solver.velocity();

                FieldSnapshot snap;
                snap.Nx = Nx;
                snap.Ny = Ny;
                snap.Lx = Lx;
                snap.Ly = Ly;
                snap.pressure.resize(static_cast<std::size_t>(Nx * Ny));
                snap.velX.resize(static_cast<std::size_t>(Nx * Ny));
                snap.velY.resize(static_cast<std::size_t>(Nx * Ny));

                for (int i = 0; i < Nx; ++i) {
                    for (int j = 0; j < Ny; ++j) {
                        const std::size_t k = static_cast<std::size_t>(i * Ny + j);
                        snap.pressure[k] = static_cast<float>(p(i, j));
                        snap.velX[k]     = static_cast<float>(u(i, j).x());
                        snap.velY[k]     = static_cast<float>(u(i, j).y());
                    }
                }
                snap.valid = true;

                std::lock_guard<std::mutex> lk(state->mtx);
                state->fieldSnap  = std::move(snap);
                state->snapReady  = true;
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

// ── letterboxView ─────────────────────────────────────────────────────────────
// Returns a view that maps the fixed logical canvas (WindowW × WindowH) into
// the actual window, maintaining aspect ratio with black bars if necessary.
static sf::View letterboxView(sf::Vector2u winSize)
{
    const float logW = UILayout::WindowW;
    const float logH = UILayout::WindowH;
    const float winW = static_cast<float>(winSize.x);
    const float winH = static_cast<float>(winSize.y);

    const float scale = std::min(winW / logW, winH / logH);
    const float vpW   = logW * scale / winW;
    const float vpH   = logH * scale / winH;
    const float vpX   = (1.f - vpW) * 0.5f;
    const float vpY   = (1.f - vpH) * 0.5f;

    sf::View view(sf::FloatRect{{0.f, 0.f}, {logW, logH}});
    view.setViewport(sf::FloatRect{{vpX, vpY}, {vpW, vpH}});
    return view;
}

// ── main ──────────────────────────────────────────────────────────────────────
int main(int /*argc*/, char* argv[])
{
    // Resolve working directory so relative paths like "cases/..." work.
    {
        namespace fs = std::filesystem;
        try {
            if (fs::exists("cases"))
                ;
#ifdef LGFLOW_SOURCE_DIR
            else if (fs::exists(fs::path(LGFLOW_SOURCE_DIR) / "cases"))
                fs::current_path(LGFLOW_SOURCE_DIR);
#endif
            else {
                fs::path dir = fs::canonical(fs::path(argv[0])).parent_path();
                for (int d = 0; d < 8; ++d) {
                    if (fs::exists(dir / "cases")) { fs::current_path(dir); break; }
                    const fs::path up = dir.parent_path();
                    if (up == dir) break;
                    dir = up;
                }
            }
        } catch (...) {}
    }

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
    window.setView(letterboxView(window.getSize()));

    // Font
    sf::Font font;
    bool fontLoaded = false;
    for (const auto& fp : std::vector<std::string>{
            "C:/Windows/Fonts/consola.ttf",
            "C:/Windows/Fonts/arial.ttf",
            "C:/Windows/Fonts/segoeui.ttf"}) {
        if (std::filesystem::exists(fp)) {
            if (font.openFromFile(fp)) { fontLoaded = true; break; }
        }
    }
    if (!fontLoaded)
        Logger::get().warn("No system font found; text may not render.");

    // UI widgets
    SidePanel       sidePanel(font);
    StatusBar       statusBar(font);
    FieldViewWidget fieldWidget(font);
    const sf::FloatRect plotBounds{
        {UILayout::PlotX, UILayout::PlotY},
        {UILayout::PlotW, UILayout::PlotH}};
    ResidualPlot plot(plotBounds, font);

    sidePanel.setCases(caseNames);
    sidePanel.setCaseDefaults({
        {32, 32,  100.0, 3000},   // Lid-Driven Cavity
        {64, 16,  100.0, 3000},   // Channel Flow
    });

    // State
    SolverState  state;
    std::thread  solverThr;
    int          selectedCase = 0;

    // Callbacks
    sidePanel.onCaseSelected = [&](int idx) {
        selectedCase = idx;
        std::lock_guard<std::mutex> lk(state.mtx);
        state.velHistory.clear();
        state.contHistory.clear();
        state.lastVel   = 0.0;
        state.lastCont  = 0.0;
        state.iteration = 0;
        state.converged = false;
        state.snapReady = false;
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
            state.snapReady = false;
        }
        plot.clear();
        const auto params = sidePanel.getParams();
        CaseOverrides ov;
        ov.Nx      = params.Nx;
        ov.Ny      = params.Ny;
        ov.Re      = params.Re;
        ov.maxIter = params.maxIter;
        const std::string cfgPath = cases[static_cast<std::size_t>(selectedCase)].second;
        if (solverThr.joinable()) solverThr.join();
        solverThr = std::thread(solverThread, cfgPath, &state, ov);
    };

    sidePanel.onStopPressed = [&]() {
        state.stopRequest.store(true);
    };

    sidePanel.onModeSelected = [&](DisplayMode mode) {
        fieldWidget.setMode(mode);
    };

    sf::Clock clock;

    while (window.isOpen()) {
        const float dt = clock.restart().asSeconds();

        while (const std::optional<sf::Event> event = window.pollEvent()) {
            if (event->is<sf::Event::Closed>()) {
                state.stopRequest.store(true);
                if (solverThr.joinable()) solverThr.join();
                window.close();
            }
            if (const auto* resized = event->getIf<sf::Event::Resized>()) {
                window.setView(letterboxView(resized->size));
            }
            sidePanel.handleEvent(*event, window);
        }

        const bool running = state.running.load();
        sidePanel.setRunning(running);

        // Update plot with new residual data
        {
            std::lock_guard<std::mutex> lk(state.mtx);
            static std::size_t lastIdx = 0;
            const std::size_t sz = state.velHistory.size();
            for (std::size_t i = lastIdx; i < sz; ++i)
                plot.push(state.velHistory[i], state.contHistory[i]);
            lastIdx = sz;
        }

        // Pull new field snapshot (outside lock for rebuild)
        {
            FieldSnapshot localSnap;
            bool hasNew = false;
            {
                std::lock_guard<std::mutex> lk(state.mtx);
                if (state.snapReady) {
                    localSnap = state.fieldSnap;
                    state.snapReady = false;
                    hasNew = true;
                }
            }
            if (hasNew)
                fieldWidget.setSnapshot(localSnap);
        }

        // Read status
        double vel = 0, cont = 0;
        int iter = 0;
        bool converged = false;
        {
            std::lock_guard<std::mutex> lk(state.mtx);
            vel       = state.lastVel;
            cont      = state.lastCont;
            iter      = state.iteration;
            converged = state.converged;
        }
        statusBar.update(iter, vel, cont, running, converged);

        const sf::Vector2f mousePos = window.mapPixelToCoords(sf::Mouse::getPosition(window));
        sidePanel.update(dt, mousePos);

        // Draw
        window.clear(sf::Color{8, 8, 12, 255});   // letterbox bars (darker than Background)
        sidePanel.draw(window);
        statusBar.draw(window);
        fieldWidget.draw(window, {{UILayout::FieldX, UILayout::FieldY},
                                  {UILayout::FieldW, UILayout::FieldH}});
        plot.draw(window);
        window.display();
    }

    state.stopRequest.store(true);
    if (solverThr.joinable()) solverThr.join();
    return 0;
}
