#pragma once
#include "RoundedRect.hpp"
#include "UIColors.hpp"
#include "UILayout.hpp"

#include <SFML/Graphics/Font.hpp>
#include <SFML/Graphics/RenderTarget.hpp>
#include <SFML/Graphics/Text.hpp>
#include <SFML/Graphics/VertexArray.hpp>

#include <algorithm>
#include <cmath>
#include <deque>
#include <string>

// ── ResidualPlot ──────────────────────────────────────────────────────────────
// Live log10 convergence plot showing vel_residual and cont_residual over
// iteration count. Draws inside the given bounds with a panel background.
class ResidualPlot
{
public:
    explicit ResidualPlot(const sf::FloatRect& bounds, const sf::Font& font)
        : m_bounds(bounds)
        , m_font(font)
    {}

    void push(double velRes, double contRes)
    {
        m_vel.push_back(static_cast<float>(velRes));
        m_cont.push_back(static_cast<float>(contRes));
        if (m_vel.size() > MAX_POINTS) {
            m_vel.pop_front();
            m_cont.pop_front();
        }
    }

    void clear()
    {
        m_vel.clear();
        m_cont.clear();
    }

    void draw(sf::RenderTarget& target) const
    {
        constexpr float pad    = UILayout::PadLg;
        const float     px     = m_bounds.position.x + pad;
        const float     py     = m_bounds.position.y + pad + 18.f;  // room for title
        const float     pw     = m_bounds.size.x - 2.f * pad;
        const float     ph     = m_bounds.size.y - 2.f * pad - 18.f;

        // Background panel
        RoundedPanel bg;
        bg.bounds       = m_bounds;
        bg.radius       = UILayout::RadiusPanel;
        bg.arcPts       = UILayout::ArcPts;
        bg.fill         = UIColors::Surface;
        bg.borderColor  = UIColors::Border;
        bg.borderThick  = UILayout::BorderThick;
        bg.shadowColor  = UIColors::Shadow;
        bg.shadowOffset = UILayout::ShadowOff;
        bg.draw(target);

        // Title
        sf::Text title(m_font, "Convergence History", 13u);
        title.setFillColor(UIColors::TextMuted);
        title.setPosition({px, m_bounds.position.y + UILayout::Pad});
        target.draw(title);

        if (m_vel.empty()) return;

        // Y axis: log10 scale from -16 to +2 (covers 1e-16 to 1e2)
        constexpr float Y_MIN = -16.f;
        constexpr float Y_MAX =   2.f;
        constexpr float Y_RNG = Y_MAX - Y_MIN;

        auto toScreenY = [&](float val) -> float {
            const float lv = (val > 0.f) ? std::log10(val) : Y_MIN;
            const float clamped = std::clamp(lv, Y_MIN, Y_MAX);
            return py + ph * (1.f - (clamped - Y_MIN) / Y_RNG);
        };
        auto toScreenX = [&](std::size_t idx, std::size_t total) -> float {
            if (total <= 1) return px;
            return px + pw * static_cast<float>(idx) / static_cast<float>(total - 1);
        };

        const std::size_t N = m_vel.size();

        // Grid lines at decades
        for (int decade = static_cast<int>(Y_MIN); decade <= static_cast<int>(Y_MAX); decade += 2) {
            const float sy = toScreenY(std::pow(10.f, static_cast<float>(decade)));
            sf::VertexArray line(sf::PrimitiveType::Lines, 2);
            line[0] = {{px,      sy}, sf::Color(60, 60, 75, 80)};
            line[1] = {{px + pw, sy}, sf::Color(60, 60, 75, 80)};
            target.draw(line);

            // Decade label
            const std::string lbl = "1e" + std::to_string(decade);
            sf::Text t(m_font, lbl, 10u);
            t.setFillColor(sf::Color(100, 100, 115, 150));
            t.setPosition({px + pw + 3.f, sy - 7.f});
            target.draw(t);
        }

        // Tolerance line at 1e-5
        {
            constexpr float tol = 1e-5f;
            const float     sy  = toScreenY(tol);
            sf::VertexArray tline(sf::PrimitiveType::Lines, 2);
            tline[0] = {{px,      sy}, sf::Color(UIColors::Success.r, UIColors::Success.g, UIColors::Success.b, 100)};
            tline[1] = {{px + pw, sy}, sf::Color(UIColors::Success.r, UIColors::Success.g, UIColors::Success.b, 100)};
            target.draw(tline);
        }

        // Plot vel_residual (accent blue)
        drawSeries(target, m_vel, N, toScreenX, toScreenY, UIColors::PlotLine);
        // Plot cont_residual (orange)
        drawSeries(target, m_cont, N, toScreenX, toScreenY, UIColors::PlotLine2);

        // Legend
        drawLegendDot(target, px,          py + ph + 4.f, UIColors::PlotLine,  "vel");
        drawLegendDot(target, px + 60.f,   py + ph + 4.f, UIColors::PlotLine2, "cont");
    }

private:
    static constexpr std::size_t MAX_POINTS = 2000;

    sf::FloatRect      m_bounds;
    const sf::Font&    m_font;
    std::deque<float>  m_vel;
    std::deque<float>  m_cont;

    template<typename FX, typename FY>
    static void drawSeries(sf::RenderTarget&          target,
                           const std::deque<float>&   data,
                           std::size_t                N,
                           FX                         toX,
                           FY                         toY,
                           sf::Color                  color)
    {
        if (N < 2) return;
        sf::VertexArray va(sf::PrimitiveType::LineStrip, N);
        for (std::size_t i = 0; i < N; ++i)
            va[i] = {{toX(i, N), toY(data[i])}, color};
        target.draw(va);
    }

    void drawLegendDot(sf::RenderTarget& target,
                       float x, float y,
                       sf::Color color,
                       const std::string& label) const
    {
        sf::VertexArray dot(sf::PrimitiveType::TriangleFan, 5);
        constexpr float r = 4.f;
        dot[0] = {{x + r, y + r}, color};
        dot[1] = {{x + 2.f * r, y + r}, color};
        dot[2] = {{x + 2.f * r, y + 2.f * r}, color};
        dot[3] = {{x, y + 2.f * r}, color};
        dot[4] = {{x, y + r}, color};
        target.draw(dot);

        sf::Text t(m_font, label, 11u);
        t.setFillColor(UIColors::TextMuted);
        t.setPosition({x + 12.f, y - 1.f});
        target.draw(t);
    }
};
