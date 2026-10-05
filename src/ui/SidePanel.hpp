#pragma once
#include "Colormap.hpp"
#include "RoundedRect.hpp"
#include "TextInput.hpp"
#include "UIColors.hpp"
#include "UILayout.hpp"

#include <SFML/Graphics/Font.hpp>
#include <SFML/Graphics/RenderTarget.hpp>
#include <SFML/Graphics/RenderWindow.hpp>
#include <SFML/Graphics/Text.hpp>
#include <SFML/Window/Event.hpp>

#include <array>
#include <cmath>
#include <functional>
#include <string>
#include <vector>

// ── AnimVar ───────────────────────────────────────────────────────────────────
// Exponential lerp variable — drives smooth hover/selection animations.
struct AnimVar
{
    float value  = 0.f;
    float target = 0.f;

    void update(float dt, float speed = 14.f)
    {
        value += (target - value) * std::min(1.f, speed * dt);
    }
};

// ── Button ────────────────────────────────────────────────────────────────────
struct Button
{
    sf::FloatRect       bounds;
    std::string         label;
    sf::Color           fillBase;
    sf::Color           fillActive;
    bool                enabled     = true;
    bool                selected    = false;

    AnimVar             hoverAnim;   // 0 = rest, 1 = hovered
    bool                hovered     = false;

    void update(float dt) { hoverAnim.update(dt); }

    void onMouseMove(sf::Vector2f pos)
    {
        hovered = enabled && bounds.contains(pos);
        hoverAnim.target = hovered ? 1.f : 0.f;
    }

    bool onMousePress(sf::Vector2f pos) const
    {
        return enabled && bounds.contains(pos);
    }

    void draw(sf::RenderTarget& target, const sf::Font& font) const
    {
        const float lift  = hoverAnim.value * 3.f;
        const float alpha = enabled ? 1.f : 0.5f;

        sf::FloatRect drawBounds{
            {bounds.position.x, bounds.position.y - lift},
            bounds.size};

        // Shadow
        RoundedPanel shadow;
        shadow.bounds       = drawBounds;
        shadow.radius       = UILayout::RadiusBtn;
        shadow.arcPts       = UILayout::ArcPts;
        shadow.fill         = UIColors::Shadow;
        shadow.shadowOffset = UILayout::ShadowOff + hoverAnim.value * 4.f;
        shadow.shadowColor  = UIColors::Shadow;
        shadow.draw(target);

        // Fill — interpolate base → active color on hover
        const auto lerp8 = [&](uint8_t a, uint8_t b) -> uint8_t {
            return static_cast<uint8_t>(a + (b - a) * hoverAnim.value);
        };
        const sf::Color col{
            lerp8(fillBase.r, fillActive.r),
            lerp8(fillBase.g, fillActive.g),
            lerp8(fillBase.b, fillActive.b),
            static_cast<uint8_t>(fillBase.a * alpha)};

        RoundedPanel btn;
        btn.bounds      = drawBounds;
        btn.radius      = UILayout::RadiusBtn;
        btn.arcPts      = UILayout::ArcPts;
        btn.fill        = col;
        if (selected) {
            btn.borderColor = UIColors::Accent;
            btn.borderThick = 1.5f;
        }
        btn.draw(target);

        // Label text — semi-transparent so it tints with button color
        sf::Text t(font, label, 13u);
        t.setFillColor(sf::Color(UIColors::TextPrimary.r,
                                 UIColors::TextPrimary.g,
                                 UIColors::TextPrimary.b, 190));
        const auto tb  = t.getLocalBounds();
        t.setPosition({
            drawBounds.position.x + (drawBounds.size.x - tb.size.x) * 0.5f - tb.position.x,
            drawBounds.position.y + (drawBounds.size.y - tb.size.y) * 0.5f - tb.position.y});
        target.draw(t);
    }
};

// ── SidePanel ─────────────────────────────────────────────────────────────────
// Left-side control panel with case selector, editable parameters, and run/stop button.
class SidePanel
{
public:
    // ── Public types ──────────────────────────────────────────────────────────
    struct CaseParams {
        int    Nx      = 32;
        int    Ny      = 32;
        double Re      = 100.0;
        int    maxIter = 3000;
    };

    // ── Callbacks ─────────────────────────────────────────────────────────────
    std::function<void(int)>         onCaseSelected;
    std::function<void()>            onRunPressed;
    std::function<void()>            onStopPressed;
    std::function<void(DisplayMode)> onModeSelected;

    explicit SidePanel(const sf::Font& font)
        : m_font(font)
    {
        rebuild();
    }

    void setCases(const std::vector<std::string>& names)
    {
        m_caseNames = names;
        rebuild();
    }

    // Call this after setCases() to supply per-case defaults.
    void setCaseDefaults(std::vector<CaseParams> defaults)
    {
        m_caseDefaults = std::move(defaults);
        if (!m_caseDefaults.empty())
            applyParams(m_caseDefaults[static_cast<std::size_t>(m_selectedCase)]);
    }

    [[nodiscard]] CaseParams getParams() const
    {
        CaseParams p;
        p.Nx      = m_inputs[0].getInt(32);
        p.Ny      = m_inputs[1].getInt(32);
        p.Re      = m_inputs[2].getDouble(100.0);
        p.maxIter = m_inputs[3].getInt(3000);
        return p;
    }

    void setRunning(bool running) { m_running = running; }

    void update(float dt, sf::Vector2f mousePos)
    {
        for (auto& btn : m_caseBtns) {
            btn.update(dt);
            btn.onMouseMove(mousePos);
        }
        for (auto& btn : m_modeBtns) {
            btn.update(dt);
            btn.onMouseMove(mousePos);
        }
        m_runBtn.update(dt);
        m_runBtn.onMouseMove(mousePos);

        for (auto& inp : m_inputs) inp.onMouseMove(mousePos);

        m_blinkTimer += dt;
        if (m_blinkTimer > 1.f) m_blinkTimer -= 1.f;
    }

    void handleEvent(const sf::Event& event, const sf::RenderWindow& window)
    {
        if (const auto* mp = event.getIf<sf::Event::MouseButtonPressed>()) {
            const sf::Vector2f pos = window.mapPixelToCoords(mp->position);
            for (int i = 0; i < static_cast<int>(m_caseBtns.size()); ++i) {
                if (m_caseBtns[static_cast<std::size_t>(i)].onMousePress(pos)) {
                    m_selectedCase = i;
                    for (auto& b : m_caseBtns) b.selected = false;
                    m_caseBtns[static_cast<std::size_t>(i)].selected = true;
                    if (!m_caseDefaults.empty() &&
                        static_cast<std::size_t>(i) < m_caseDefaults.size())
                        applyParams(m_caseDefaults[static_cast<std::size_t>(i)]);
                    if (onCaseSelected) onCaseSelected(i);
                }
            }
            if (m_runBtn.onMousePress(pos)) {
                if (m_running) { if (onStopPressed) onStopPressed(); }
                else           { if (onRunPressed)  onRunPressed();  }
            }
            const DisplayMode modes[] = {
                DisplayMode::Pressure, DisplayMode::Velocity, DisplayMode::Streamlines};
            for (int i = 0; i < static_cast<int>(m_modeBtns.size()); ++i) {
                if (m_modeBtns[static_cast<std::size_t>(i)].onMousePress(pos)) {
                    m_displayMode = modes[i];
                    for (auto& b : m_modeBtns) b.selected = false;
                    m_modeBtns[static_cast<std::size_t>(i)].selected = true;
                    if (onModeSelected) onModeSelected(m_displayMode);
                }
            }
            // TextInput focus handling — blur all then let the clicked one focus
            bool anyFocused = false;
            for (auto& inp : m_inputs) {
                if (inp.onMousePress(pos)) anyFocused = true;
            }
            // If a different input was clicked, the others are already blurred
            // because onMousePress() sets focused = contains(pos).
            (void)anyFocused;
        }
        if (const auto* mm = event.getIf<sf::Event::MouseMoved>()) {
            const sf::Vector2f pos = window.mapPixelToCoords(mm->position);
            for (auto& btn : m_caseBtns) btn.onMouseMove(pos);
            for (auto& btn : m_modeBtns) btn.onMouseMove(pos);
            m_runBtn.onMouseMove(pos);
        }
        if (const auto* te = event.getIf<sf::Event::TextEntered>()) {
            for (auto& inp : m_inputs) inp.handleText(te->unicode);
        }
        if (const auto* kp = event.getIf<sf::Event::KeyPressed>()) {
            for (auto& inp : m_inputs) inp.handleKey(kp->code);
        }
    }

    void draw(sf::RenderTarget& target) const
    {
        // Panel background
        RoundedPanel bg;
        bg.bounds       = {
            {UILayout::SideX, UILayout::SideY},
            {UILayout::SideW, UILayout::SideH}};
        bg.radius       = UILayout::RadiusPanel;
        bg.arcPts       = UILayout::ArcPts;
        bg.fill         = UIColors::Surface;
        bg.borderColor  = UIColors::Border;
        bg.borderThick  = UILayout::BorderThick;
        bg.shadowColor  = UIColors::Shadow;
        bg.shadowOffset = UILayout::ShadowOff;
        bg.draw(target);

        // Title
        sf::Text title(m_font, "LG-Flow", 18u);
        title.setFillColor(UIColors::Accent);
        title.setPosition({UILayout::SideX + UILayout::Pad,
                           UILayout::SideY + UILayout::Pad});
        target.draw(title);

        sf::Text sub(m_font, "CFD Solver", 11u);
        sub.setFillColor(UIColors::TextMuted);
        sub.setPosition({UILayout::SideX + UILayout::Pad,
                         UILayout::SideY + UILayout::Pad + 22.f});
        target.draw(sub);

        // Section label — cases
        sf::Text lbl(m_font, "CASES", 10u);
        lbl.setFillColor(UIColors::TextMuted);
        lbl.setPosition({UILayout::SideX + UILayout::Pad,
                         UILayout::SideY + 66.f});
        target.draw(lbl);

        // Case buttons
        for (const auto& btn : m_caseBtns)
            btn.draw(target, m_font);

        // Section label — params
        sf::Text plbl(m_font, "PARAMS", 10u);
        plbl.setFillColor(UIColors::TextMuted);
        plbl.setPosition({UILayout::SideX + UILayout::Pad, m_paramsLabelY});
        target.draw(plbl);

        // Parameter inputs
        for (const auto& inp : m_inputs)
            inp.draw(target, m_font, m_blinkTimer);

        // Section label — display mode
        sf::Text dlbl(m_font, "DISPLAY", 10u);
        dlbl.setFillColor(UIColors::TextMuted);
        dlbl.setPosition({UILayout::SideX + UILayout::Pad, m_displayLabelY});
        target.draw(dlbl);

        // Display mode buttons
        for (const auto& btn : m_modeBtns)
            btn.draw(target, m_font);

        // Run/Stop button
        Button runDraw = m_runBtn;
        runDraw.label     = m_running ? "Stop" : "Run";
        runDraw.fillBase  = m_running ? sf::Color{80, 40, 40, 255}
                                      : sf::Color{40, 80, 50, 255};
        runDraw.fillActive = m_running ? UIColors::Danger : UIColors::Success;
        runDraw.draw(target, m_font);
    }

private:
    const sf::Font&           m_font;
    std::vector<std::string>  m_caseNames;
    std::vector<Button>       m_caseBtns;
    std::vector<Button>       m_modeBtns;
    Button                    m_runBtn;
    std::array<TextInput, 4>  m_inputs;
    std::vector<CaseParams>   m_caseDefaults;
    int                       m_selectedCase  = 0;
    bool                      m_running       = false;
    DisplayMode               m_displayMode   = DisplayMode::Pressure;
    float                     m_displayLabelY = 0.f;
    float                     m_paramsLabelY  = 0.f;
    float                     m_blinkTimer    = 0.f;

    void applyParams(const CaseParams& p)
    {
        m_inputs[0].setValue(p.Nx);
        m_inputs[1].setValue(p.Ny);
        m_inputs[2].setValue(p.Re);
        m_inputs[3].setValue(p.maxIter);
    }

    void rebuild()
    {
        m_caseBtns.clear();
        m_modeBtns.clear();

        const float bx = UILayout::SideX + UILayout::Pad;
        const float bw = UILayout::SideW - 2.f * UILayout::Pad;
        constexpr float ModeH   = 34.f;
        constexpr float InputH  = 26.f;
        constexpr float InputGap =  5.f;

        // Case buttons
        float y = UILayout::SideY + 84.f;
        for (std::size_t i = 0; i < m_caseNames.size(); ++i) {
            Button btn;
            btn.bounds     = {{bx, y}, {bw, UILayout::BtnH}};
            btn.label      = m_caseNames[i];
            btn.fillBase   = UIColors::SurfaceHigh;
            btn.fillActive = UIColors::AccentDim;
            btn.selected   = (static_cast<int>(i) == m_selectedCase);
            m_caseBtns.push_back(btn);
            y += UILayout::BtnH + UILayout::Gap;
        }

        // PARAMS section
        m_paramsLabelY = y + UILayout::Gap;
        y = m_paramsLabelY + 16.f;

        const char* inputLabels[] = {"Nx", "Ny", "Re", "Max Iter"};
        const bool  inputIntOnly[] = {true, true, false, true};
        for (int i = 0; i < 4; ++i) {
            m_inputs[static_cast<std::size_t>(i)].bounds  = {{bx, y}, {bw, InputH}};
            m_inputs[static_cast<std::size_t>(i)].label   = inputLabels[i];
            m_inputs[static_cast<std::size_t>(i)].intOnly = inputIntOnly[i];
            y += InputH + InputGap;
        }

        // Apply defaults for the selected case if available
        if (!m_caseDefaults.empty() &&
            static_cast<std::size_t>(m_selectedCase) < m_caseDefaults.size())
            applyParams(m_caseDefaults[static_cast<std::size_t>(m_selectedCase)]);

        // Display mode buttons
        m_displayLabelY = y + UILayout::Gap;
        y = m_displayLabelY + 16.f;
        const char* modeLabels[] = {"Pressure", "Velocity", "Streamlines"};
        const DisplayMode modes[] = {
            DisplayMode::Pressure, DisplayMode::Velocity, DisplayMode::Streamlines};
        for (int i = 0; i < 3; ++i) {
            Button btn;
            btn.bounds     = {{bx, y}, {bw, ModeH}};
            btn.label      = modeLabels[i];
            btn.fillBase   = UIColors::SurfaceHigh;
            btn.fillActive = UIColors::AccentDim;
            btn.selected   = (modes[i] == m_displayMode);
            m_modeBtns.push_back(btn);
            y += ModeH + UILayout::Gap;
        }

        // Run button at bottom
        const float runY = UILayout::SideY + UILayout::SideH
                         - UILayout::BtnH - UILayout::Pad;
        m_runBtn.bounds    = {{bx, runY}, {bw, UILayout::BtnH}};
        m_runBtn.label     = "Run";
        m_runBtn.fillBase  = sf::Color{40, 80, 50, 255};
        m_runBtn.fillActive = UIColors::Success;
    }
};
