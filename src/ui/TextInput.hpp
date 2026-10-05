#pragma once
#include "RoundedRect.hpp"
#include "UIColors.hpp"
#include "UILayout.hpp"

#include <SFML/Graphics/Font.hpp>
#include <SFML/Graphics/RenderTarget.hpp>
#include <SFML/Graphics/Text.hpp>
#include <SFML/Graphics/VertexArray.hpp>
#include <SFML/Window/Event.hpp>
#include <SFML/Window/Keyboard.hpp>

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <string>

// ── TextInput ─────────────────────────────────────────────────────────────────
// Single-line text entry field with label, digit-only filter, and focus state.
struct TextInput
{
    sf::FloatRect bounds;
    std::string   label;     ///< Short label drawn left of the input box
    std::string   value;     ///< Current text contents
    bool          intOnly = true;  ///< Accept only digits (no '.' or 'e')
    bool          focused  = false;
    bool          hovered  = false;

    [[nodiscard]] int    getInt(int fallback = 0) const
    {
        try { return std::stoi(value); } catch (...) { return fallback; }
    }
    [[nodiscard]] double getDouble(double fallback = 0.0) const
    {
        try { return std::stod(value); } catch (...) { return fallback; }
    }

    void setValue(int v)    { value = std::to_string(v); }
    void setValue(double v)
    {
        std::ostringstream ss;
        ss << v;
        value = ss.str();
    }

    void onMouseMove(sf::Vector2f pos)
    {
        hovered = bounds.contains(pos);
    }

    // Returns true if this field claimed focus.
    bool onMousePress(sf::Vector2f pos)
    {
        focused = bounds.contains(pos);
        return focused;
    }

    void blur() { focused = false; }

    void handleText(uint32_t unicode)
    {
        if (!focused) return;
        if (unicode == 8u) {  // Backspace
            if (!value.empty()) value.pop_back();
            return;
        }
        if (unicode < 32u || unicode > 126u) return;
        const char c = static_cast<char>(unicode);
        if (intOnly) {
            if (c >= '0' && c <= '9') value += c;
        } else {
            if ((c >= '0' && c <= '9') || c == '.' || c == '-' || c == 'e')
                value += c;
        }
    }

    void handleKey(sf::Keyboard::Key key)
    {
        if (!focused) return;
        if (key == sf::Keyboard::Key::Enter || key == sf::Keyboard::Key::Escape)
            focused = false;
    }

    void draw(sf::RenderTarget& target, const sf::Font& font,
              float blinkPhase) const
    {
        const float LabelW = 52.f;   // width reserved for the label text
        const float gap    =  6.f;

        // Label
        sf::Text lbl(font, label, 11u);
        lbl.setFillColor(UIColors::TextMuted);
        const auto lb = lbl.getLocalBounds();
        lbl.setPosition({
            bounds.position.x + LabelW - lb.size.x - lb.position.x,
            bounds.position.y + (bounds.size.y - lb.size.y) * 0.5f - lb.position.y});
        target.draw(lbl);

        // Input box background
        const sf::FloatRect boxBounds{
            {bounds.position.x + LabelW + gap, bounds.position.y},
            {bounds.size.x - LabelW - gap, bounds.size.y}};

        RoundedPanel box;
        box.bounds      = boxBounds;
        box.radius      = UILayout::RadiusInner;
        box.arcPts      = UILayout::ArcPts;
        box.fill        = focused  ? UIColors::AccentDim
                        : hovered  ? UIColors::SurfaceHigh
                                   : sf::Color{28, 28, 38, 255};
        box.borderColor = focused  ? UIColors::Accent : UIColors::Border;
        box.borderThick = focused  ? 1.5f : 1.f;
        box.draw(target);

        // Value text (right-aligned inside box, with space for cursor)
        constexpr float TextPad = 5.f;
        sf::Text val(font, value, 11u);
        val.setFillColor(UIColors::TextPrimary);
        const auto vb = val.getLocalBounds();
        val.setPosition({
            boxBounds.position.x + boxBounds.size.x - vb.size.x - vb.position.x - TextPad - (focused ? 6.f : 0.f),
            boxBounds.position.y + (boxBounds.size.y - vb.size.y) * 0.5f - vb.position.y});
        target.draw(val);

        // Cursor (blinking, shown only when focused)
        if (focused && blinkPhase < 0.5f) {
            const float cx = boxBounds.position.x + boxBounds.size.x - TextPad - 4.f;
            const float cy = boxBounds.position.y + 4.f;
            const float ch = boxBounds.size.y - 8.f;
            sf::VertexArray cursor(sf::PrimitiveType::Lines, 2);
            cursor[0] = {{cx, cy},      UIColors::Accent};
            cursor[1] = {{cx, cy + ch}, UIColors::Accent};
            target.draw(cursor);
        }
    }
};
