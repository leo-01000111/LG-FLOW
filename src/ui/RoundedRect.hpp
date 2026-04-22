#pragma once
#include <SFML/Graphics/Color.hpp>
#include <SFML/Graphics/RenderTarget.hpp>
#include <SFML/Graphics/Vertex.hpp>
#include <SFML/Graphics/VertexArray.hpp>
#include <SFML/System/Vector2.hpp>

#include <cmath>
#include <numbers>
#include <vector>

// ── CircleArcGenerator ────────────────────────────────────────────────────────
// Generates arc points for one quarter-circle, used internally by
// RoundedRectGenerator. pointCount includes both endpoints.
struct CircleArcGenerator
{
    sf::Vector2f center;
    float        radius;
    float        startAngleDeg;
    int          pointCount;

    [[nodiscard]] std::vector<sf::Vector2f> generate() const
    {
        std::vector<sf::Vector2f> pts;
        pts.reserve(static_cast<std::size_t>(pointCount));
        for (int i = 0; i < pointCount; ++i) {
            const float t   = static_cast<float>(i) / static_cast<float>(pointCount - 1);
            const float deg = startAngleDeg + t * 90.f;
            const float rad = deg * (std::numbers::pi_v<float> / 180.f);
            pts.push_back(center + sf::Vector2f(std::cos(rad), std::sin(rad)) * radius);
        }
        return pts;
    }
};

// ── RoundedRectGenerator ──────────────────────────────────────────────────────
// Produces the outline points of a rounded rectangle as a closed polygon.
// Uses 4 × CircleArcGenerator; seam points are not duplicated.
struct RoundedRectGenerator
{
    sf::FloatRect bounds;
    float         cornerRadius;
    int           pointsPerCorner;  // arc quality knob

    [[nodiscard]] std::vector<sf::Vector2f> generate() const
    {
        const float r  = cornerRadius;
        const float l  = bounds.position.x;
        const float t  = bounds.position.y;
        const float ri = l + bounds.size.x;
        const float bo = t + bounds.size.y;

        // Corner centres
        const sf::Vector2f tl{l + r, t + r};
        const sf::Vector2f tr{ri - r, t + r};
        const sf::Vector2f br{ri - r, bo - r};
        const sf::Vector2f bl{l + r, bo - r};

        // Each arc uses (pointsPerCorner - 1) arcs so shared corner points
        // between adjacent arcs are not duplicated → total = 4*(pts-1).
        const int pts = pointsPerCorner;
        std::vector<sf::Vector2f> poly;
        poly.reserve(static_cast<std::size_t>(4 * (pts - 1)));

        auto appendArc = [&](sf::Vector2f c, float startDeg) {
            CircleArcGenerator arc{c, r, startDeg, pts};
            auto ap = arc.generate();
            // Drop the last point — it is shared with the next arc's first point
            for (int i = 0; i < pts - 1; ++i)
                poly.push_back(ap[static_cast<std::size_t>(i)]);
        };

        appendArc(tl, 180.f);  // top-left:     180 → 270
        appendArc(tr, 270.f);  // top-right:    270 → 360
        appendArc(br,   0.f);  // bottom-right:   0 →  90
        appendArc(bl,  90.f);  // bottom-left:   90 → 180

        return poly;
    }
};

// ── Helper: build a filled TriangleFan VertexArray from outline points ─────────
inline sf::VertexArray makeFilled(const std::vector<sf::Vector2f>& outline,
                                  sf::Color                         fill)
{
    const auto N = static_cast<std::size_t>(outline.size());
    sf::VertexArray va(sf::PrimitiveType::TriangleFan, N + 2);

    // Hub = centroid (average of all points)
    sf::Vector2f hub{0.f, 0.f};
    for (const auto& p : outline) hub += p;
    hub /= static_cast<float>(N);

    va[0] = {hub, fill};
    for (std::size_t i = 0; i < N; ++i)
        va[i + 1] = {outline[i], fill};
    va[N + 1] = {outline[0], fill};  // close the fan

    return va;
}

// ── Helper: build an outline TriangleStrip (double-generator) ─────────────────
inline sf::VertexArray makeOutline(const std::vector<sf::Vector2f>& outer,
                                   const std::vector<sf::Vector2f>& inner,
                                   sf::Color                         color)
{
    const auto N = outer.size();
    sf::VertexArray va(sf::PrimitiveType::TriangleStrip, 2 * N + 2);
    for (std::size_t i = 0; i < N; ++i) {
        va[2 * i]     = {outer[i], color};
        va[2 * i + 1] = {inner[i], color};
    }
    // Close the strip
    va[2 * N]     = {outer[0], color};
    va[2 * N + 1] = {inner[0], color};
    return va;
}

// ── RoundedPanel ─────────────────────────────────────────────────────────────
// Convenience: draws a filled rounded rect + optional drop shadow + optional border.
// Follows the style guide's layout build order (depth before color).
struct RoundedPanel
{
    sf::FloatRect bounds;
    float         radius;
    int           arcPts;
    sf::Color     fill;
    sf::Color     borderColor  = sf::Color::Transparent;
    float         borderThick  = 0.f;
    sf::Color     shadowColor  = sf::Color::Transparent;
    float         shadowOffset = 0.f;

    void draw(sf::RenderTarget& target) const
    {
        RoundedRectGenerator gen{bounds, radius, arcPts};

        // Shadow (offset downward)
        if (shadowOffset > 0.f && shadowColor.a > 0) {
            const sf::FloatRect sb{
                {bounds.position.x, bounds.position.y + shadowOffset},
                bounds.size};
            RoundedRectGenerator sg{sb, radius, arcPts};
            auto sv = makeFilled(sg.generate(), shadowColor);
            target.draw(sv);
        }

        // Fill
        auto fv = makeFilled(gen.generate(), fill);
        target.draw(fv);

        // Border
        if (borderThick > 0.f && borderColor.a > 0) {
            const sf::FloatRect inner{
                {bounds.position.x + borderThick,
                 bounds.position.y + borderThick},
                {bounds.size.x - 2.f * borderThick,
                 bounds.size.y - 2.f * borderThick}};
            const float innerR = std::max(0.f, radius - borderThick);
            RoundedRectGenerator ig{inner, innerR, arcPts};
            auto bv = makeOutline(gen.generate(), ig.generate(), borderColor);
            target.draw(bv);
        }
    }
};
