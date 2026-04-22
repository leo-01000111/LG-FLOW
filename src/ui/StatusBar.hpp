#pragma once
#include "RoundedRect.hpp"
#include "UIColors.hpp"
#include "UILayout.hpp"

#include <SFML/Graphics/Font.hpp>
#include <SFML/Graphics/RenderTarget.hpp>
#include <SFML/Graphics/Text.hpp>

#include <format>
#include <string>

// ── StatusBar ─────────────────────────────────────────────────────────────────
// Thin bar at the top-right showing iteration count, residuals, and state pill.
class StatusBar
{
public:
    explicit StatusBar(const sf::Font& font) : m_font(font) {}

    void update(int iter, double velRes, double contRes, bool running, bool converged)
    {
        m_iter      = iter;
        m_vel       = velRes;
        m_cont      = contRes;
        m_running   = running;
        m_converged = converged;
    }

    void draw(sf::RenderTarget& target) const
    {
        const sf::FloatRect bounds{
            {UILayout::StatusX, UILayout::StatusY},
            {UILayout::StatusW, UILayout::StatusH}};

        RoundedPanel bg;
        bg.bounds       = bounds;
        bg.radius       = UILayout::RadiusCard;
        bg.arcPts       = UILayout::ArcPts;
        bg.fill         = UIColors::Surface;
        bg.borderColor  = UIColors::Border;
        bg.borderThick  = UILayout::BorderThick;
        bg.shadowColor  = UIColors::Shadow;
        bg.shadowOffset = UILayout::ShadowOff;
        bg.draw(target);

        const float cx = bounds.position.x + UILayout::Pad;
        const float cy = bounds.position.y + (bounds.size.y - 14.f) * 0.5f;

        // State pill
        const sf::Color pillCol = m_converged ? UIColors::Success
                                 : m_running  ? UIColors::Accent
                                              : UIColors::Border;
        const std::string pillTxt = m_converged ? "Converged"
                                   : m_running   ? "Running"
                                                 : "Ready";
        drawPill(target, cx, cy, pillTxt, pillCol);

        // Stats
        const float statsX = cx + 100.f;
        drawStat(target, statsX,        cy, "iter", std::to_string(m_iter));
        drawStat(target, statsX + 110.f, cy, "vel",  fmtSci(m_vel));
        drawStat(target, statsX + 260.f, cy, "cont", fmtSci(m_cont));
    }

private:
    const sf::Font& m_font;
    int    m_iter      = 0;
    double m_vel       = 0.0;
    double m_cont      = 0.0;
    bool   m_running   = false;
    bool   m_converged = false;

    static std::string fmtSci(double v)
    {
        if (v == 0.0) return "0";
        return std::format("{:.2e}", v);
    }

    void drawPill(sf::RenderTarget& target,
                  float x, float y,
                  const std::string& text,
                  sf::Color color) const
    {
        sf::Text t(m_font, text, 12u);
        const auto tb = t.getLocalBounds();
        const float pw = tb.size.x + 16.f;
        const float ph = 22.f;
        const float pr = ph * 0.5f;

        RoundedPanel pill;
        pill.bounds   = {{x, y}, {pw, ph}};
        pill.radius   = pr;
        pill.arcPts   = UILayout::ArcPts;
        pill.fill     = sf::Color(color.r, color.g, color.b, 40);
        pill.borderColor = color;
        pill.borderThick = 1.f;
        pill.draw(target);

        t.setFillColor(color);
        t.setPosition({x + 8.f - tb.position.x,
                       y + (ph - tb.size.y) * 0.5f - tb.position.y});
        target.draw(t);
    }

    void drawStat(sf::RenderTarget& target,
                  float x, float y,
                  const std::string& key,
                  const std::string& val) const
    {
        sf::Text k(m_font, key + ":", 11u);
        k.setFillColor(UIColors::TextMuted);
        k.setPosition({x, y + 2.f});
        target.draw(k);

        const auto kb = k.getLocalBounds();
        sf::Text v(m_font, val, 12u);
        v.setFillColor(UIColors::TextPrimary);
        v.setPosition({x + kb.size.x + kb.position.x + 4.f, y + 1.f});
        target.draw(v);
    }
};
