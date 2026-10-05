#pragma once
#include "Colormap.hpp"
#include "RoundedRect.hpp"
#include "UIColors.hpp"
#include "UILayout.hpp"

#include <SFML/Graphics/Font.hpp>
#include <SFML/Graphics/RenderTarget.hpp>
#include <SFML/Graphics/Text.hpp>
#include <SFML/Graphics/VertexArray.hpp>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

// ── FieldSnapshot ─────────────────────────────────────────────────────────────
struct FieldSnapshot
{
    std::vector<float> pressure;  // [i*Ny + j]
    std::vector<float> velX;
    std::vector<float> velY;
    int   Nx = 0, Ny = 0;
    float Lx = 0.f, Ly = 0.f;
    bool  valid = false;
};

// ── FieldView ─────────────────────────────────────────────────────────────────
// Renders pressure/velocity heatmaps and streamlines inside a given bounds rect.
class FieldView
{
public:
    explicit FieldView(const sf::Font& font) : m_font(font) {}

    void setMode(DisplayMode mode)
    {
        m_mode = mode;
        if (m_snap.valid) rebuild();
    }
    DisplayMode mode() const { return m_mode; }

    void setSnapshot(const FieldSnapshot& snap)
    {
        m_snap = snap;
        if (m_snap.valid) rebuild();
    }

    void draw(sf::RenderTarget& target, sf::FloatRect bounds) const
    {
        // Panel background
        RoundedPanel bg;
        bg.bounds       = bounds;
        bg.radius       = UILayout::RadiusPanel;
        bg.arcPts       = UILayout::ArcPts;
        bg.fill         = UIColors::Surface;
        bg.borderColor  = UIColors::Border;
        bg.borderThick  = UILayout::BorderThick;
        bg.shadowColor  = UIColors::Shadow;
        bg.shadowOffset = UILayout::ShadowOff;
        bg.draw(target);

        // Title
        const std::string titles[] = {"Pressure", "Velocity Magnitude", "Streamlines"};
        sf::Text title(m_font, titles[static_cast<int>(m_mode)], 13u);
        title.setFillColor(UIColors::TextMuted);
        title.setPosition({bounds.position.x + UILayout::Pad,
                           bounds.position.y + UILayout::Pad});
        target.draw(title);

        if (!m_snap.valid) {
            sf::Text hint(m_font, "Run a case to see the flow field", 13u);
            hint.setFillColor(UIColors::TextMuted);
            const auto hb = hint.getLocalBounds();
            hint.setPosition({
                bounds.position.x + (bounds.size.x - hb.size.x) * 0.5f,
                bounds.position.y + bounds.size.y * 0.5f});
            target.draw(hint);
            return;
        }

        // Inner draw area: leave room for title (top) and colorbar (right)
        constexpr float TitleH  = 24.f;
        constexpr float CbarW   = 28.f;
        constexpr float CbarGap = 6.f;
        const sf::FloatRect inner{
            {bounds.position.x + UILayout::Pad,
             bounds.position.y + UILayout::Pad + TitleH},
            {bounds.size.x - 2.f * UILayout::Pad - CbarW - CbarGap,
             bounds.size.y - 2.f * UILayout::Pad - TitleH}};

        // Coordinate transform: physical → screen (with y-flip)
        const float sx = inner.position.x;
        const float sy = inner.position.y;
        const float sw = inner.size.x;
        const float sh = inner.size.y;
        const float Lx = m_snap.Lx;
        const float Ly = m_snap.Ly;

        auto toScr = [&](float px, float py) -> sf::Vector2f {
            return {sx + px / Lx * sw,
                    sy + sh - py / Ly * sh};
        };

        // Draw heatmap grid
        target.draw(m_heatmap);

        // Draw streamlines (mode == Streamlines)
        if (m_mode == DisplayMode::Streamlines) {
            for (const auto& sl : m_streamlines) {
                if (sl.size() < 2) continue;
                sf::VertexArray va(sf::PrimitiveType::LineStrip, sl.size());
                for (std::size_t i = 0; i < sl.size(); ++i) {
                    const sf::Vector2f s = toScr(sl[i].x, sl[i].y);
                    va[i] = {s, sl[i].color};
                }
                target.draw(va);
            }
        }

        // Draw velocity arrows (pressure and velocity modes)
        if (m_mode != DisplayMode::Streamlines) {
            target.draw(m_arrows);
        }

        // Aspect-ratio border (thin rect outline around field area)
        {
            const sf::Vector2f bl = toScr(0.f,  0.f);
            const sf::Vector2f tr = toScr(Lx, Ly);
            sf::VertexArray border(sf::PrimitiveType::LineStrip, 5);
            border[0] = {{tr.x, sy}, UIColors::Border};
            border[1] = {{bl.x, sy}, UIColors::Border};
            border[2] = {{bl.x, bl.y}, UIColors::Border};
            border[3] = {{tr.x, bl.y}, UIColors::Border};
            border[4] = {{tr.x, sy}, UIColors::Border};
            target.draw(border);
        }

        // Colorbar
        {
            const float cbX = inner.position.x + inner.size.x + CbarGap;
            const float cbY = inner.position.y;
            const float cbH = inner.size.y;
            constexpr int N = 64;
            sf::VertexArray cb(sf::PrimitiveType::TriangleStrip, static_cast<std::size_t>(2 * N));
            for (int k = 0; k < N; ++k) {
                const float t  = static_cast<float>(k) / static_cast<float>(N - 1);
                const float y0 = cbY + cbH * (1.f - t);
                const sf::Color col = (m_mode == DisplayMode::Velocity)
                                    ? Colormap::velocity(t)
                                    : Colormap::pressure(t);
                cb[static_cast<std::size_t>(2 * k)]     = {{cbX,           y0}, col};
                cb[static_cast<std::size_t>(2 * k + 1)] = {{cbX + CbarW,   y0}, col};
            }
            target.draw(cb);

            // Min/max labels
            sf::Text tmax(m_font, fmtVal(m_scalarMax), 10u);
            tmax.setFillColor(UIColors::TextMuted);
            tmax.setPosition({cbX, cbY - 1.f});
            target.draw(tmax);

            sf::Text tmin(m_font, fmtVal(m_scalarMin), 10u);
            tmin.setFillColor(UIColors::TextMuted);
            const auto tb = tmin.getLocalBounds();
            tmin.setPosition({cbX, cbY + cbH - tb.size.y - tb.position.y});
            target.draw(tmin);
        }
    }

private:
    // A streamline point carries its color (from speed at that location)
    struct SLPoint { float x, y; sf::Color color; };

    const sf::Font&             m_font;
    DisplayMode                 m_mode = DisplayMode::Pressure;
    FieldSnapshot               m_snap;

    sf::VertexArray             m_heatmap;   // 6 vertices per cell (2 triangles)
    sf::VertexArray             m_arrows;    // velocity arrow lines
    std::vector<std::vector<SLPoint>> m_streamlines;

    float m_scalarMin = 0.f;
    float m_scalarMax = 1.f;

    static std::string fmtVal(float v)
    {
        char buf[32];
        if (std::abs(v) < 1e-3f || std::abs(v) > 999.f)
            std::snprintf(buf, sizeof(buf), "%.1e", static_cast<double>(v));
        else
            std::snprintf(buf, sizeof(buf), "%.3f", static_cast<double>(v));
        return buf;
    }

    // Bilinear interpolation of velocity at physical position (px, py)
    sf::Vector2f interpVel(float px, float py) const
    {
        const float dx = m_snap.Lx / static_cast<float>(m_snap.Nx);
        const float dy = m_snap.Ly / static_cast<float>(m_snap.Ny);
        const int Nx = m_snap.Nx, Ny = m_snap.Ny;

        float fi = px / dx - 0.5f;
        float fj = py / dy - 0.5f;
        fi = std::clamp(fi, 0.f, static_cast<float>(Nx - 1));
        fj = std::clamp(fj, 0.f, static_cast<float>(Ny - 1));

        const int i0 = std::clamp(static_cast<int>(fi),     0, Nx - 1);
        const int i1 = std::clamp(static_cast<int>(fi) + 1, 0, Nx - 1);
        const int j0 = std::clamp(static_cast<int>(fj),     0, Ny - 1);
        const int j1 = std::clamp(static_cast<int>(fj) + 1, 0, Ny - 1);
        const float tx = fi - static_cast<float>(i0);
        const float ty = fj - static_cast<float>(j0);

        auto idx = [&](int i, int j) { return static_cast<std::size_t>(i * Ny + j); };
        const float w00 = (1.f - tx) * (1.f - ty);
        const float w10 = tx * (1.f - ty);
        const float w01 = (1.f - tx) * ty;
        const float w11 = tx * ty;

        return {
            w00 * m_snap.velX[idx(i0,j0)] + w10 * m_snap.velX[idx(i1,j0)]
          + w01 * m_snap.velX[idx(i0,j1)] + w11 * m_snap.velX[idx(i1,j1)],
            w00 * m_snap.velY[idx(i0,j0)] + w10 * m_snap.velY[idx(i1,j0)]
          + w01 * m_snap.velY[idx(i0,j1)] + w11 * m_snap.velY[idx(i1,j1)]};
    }

    // Integrate one streamline using arc-length-parametrised RK4.
    // Returns list of SLPoint (physical coords + speed color).
    std::vector<SLPoint> integrateStreamline(float x0, float y0,
                                             float velMax) const
    {
        const float ds = 0.15f * std::min(m_snap.Lx / m_snap.Nx,
                                          m_snap.Ly / m_snap.Ny);
        constexpr int MAX_STEPS = 2000;

        std::vector<SLPoint> pts;
        pts.reserve(256);
        float x = x0, y = y0;

        for (int step = 0; step < MAX_STEPS; ++step) {
            // Check bounds
            if (x < 0.f || x > m_snap.Lx || y < 0.f || y > m_snap.Ly) break;

            const sf::Vector2f v = interpVel(x, y);
            const float spd = std::sqrt(v.x * v.x + v.y * v.y);
            if (spd < 1e-6f) break;

            const float t = (velMax > 1e-6f) ? spd / velMax : 0.f;
            pts.push_back({x, y, Colormap::velocity(t)});

            // Normalised direction for arc-length step (RK4)
            auto dir = [&](float px, float py) -> sf::Vector2f {
                const sf::Vector2f vv = interpVel(px, py);
                const float s = std::sqrt(vv.x * vv.x + vv.y * vv.y);
                if (s < 1e-9f) return {0.f, 0.f};
                return {vv.x / s, vv.y / s};
            };

            const sf::Vector2f k1 = dir(x,                  y);
            const sf::Vector2f k2 = dir(x + 0.5f*ds*k1.x,  y + 0.5f*ds*k1.y);
            const sf::Vector2f k3 = dir(x + 0.5f*ds*k2.x,  y + 0.5f*ds*k2.y);
            const sf::Vector2f k4 = dir(x +      ds*k3.x,  y +      ds*k3.y);

            x += ds / 6.f * (k1.x + 2.f*k2.x + 2.f*k3.x + k4.x);
            y += ds / 6.f * (k1.y + 2.f*k2.y + 2.f*k3.y + k4.y);
        }
        return pts;
    }

    void rebuild()
    {
        rebuildHeatmap();
        rebuildArrows();
        if (m_mode == DisplayMode::Streamlines) rebuildStreamlines();
    }

    void rebuildHeatmap()
    {
        const int Nx = m_snap.Nx, Ny = m_snap.Ny;
        const int Nc = Nx * Ny;

        // Compute scalar field and range
        std::vector<float> scalar(static_cast<std::size_t>(Nc));
        if (m_mode == DisplayMode::Velocity || m_mode == DisplayMode::Streamlines) {
            for (int k = 0; k < Nc; ++k)
                scalar[static_cast<std::size_t>(k)] = std::sqrt(
                    m_snap.velX[static_cast<std::size_t>(k)] * m_snap.velX[static_cast<std::size_t>(k)] +
                    m_snap.velY[static_cast<std::size_t>(k)] * m_snap.velY[static_cast<std::size_t>(k)]);
        } else {
            scalar = m_snap.pressure;
        }

        m_scalarMin = *std::min_element(scalar.begin(), scalar.end());
        m_scalarMax = *std::max_element(scalar.begin(), scalar.end());
        if (m_scalarMax - m_scalarMin < 1e-10f) m_scalarMax = m_scalarMin + 1e-10f;

        // Build vertex array in NORMALISED coordinates [0,1]x[0,1].
        // Actual screen positions are set in draw() via a uniform scale,
        // but we embed physical fractions directly as pixel positions scaled
        // by unit bounds — we'll let SFML draw with a view transform.
        //
        // Simpler: store actual physical fractions, scale in draw() via
        // the sf::Transform approach.  Here we store UNIT coordinates
        // in [0, Nx] x [0, Ny] and map in draw() with a view.
        //
        // ACTUALLY — for simplicity and to avoid view complications,
        // store vertices in "fraction of field" units [0,1]x[0,1].
        // draw() will apply the affine transform to screen space.
        // But VertexArray doesn't support per-draw transforms directly.
        //
        // FINAL decision: store in absolute screen pixels with a dummy
        // bounds [0,1] and rescale in draw() ... no.
        //
        // Cleanest for a header-only class: store in physical units
        // and draw with sf::RenderStates + Transform.
        // Compute the transform in draw() and apply here as a pre-transform.
        //
        // Here we store PHYSICAL coords; draw() computes the transform.
        m_heatmap = sf::VertexArray(sf::PrimitiveType::Triangles,
                                    static_cast<std::size_t>(Nc * 6));

        const float dx = m_snap.Lx / static_cast<float>(Nx);
        const float dy = m_snap.Ly / static_cast<float>(Ny);

        for (int i = 0; i < Nx; ++i) {
            for (int j = 0; j < Ny; ++j) {
                const float t = (scalar[static_cast<std::size_t>(i * Ny + j)] - m_scalarMin)
                              / (m_scalarMax - m_scalarMin);
                const sf::Color col = (m_mode == DisplayMode::Pressure)
                                    ? Colormap::pressure(t)
                                    : Colormap::velocity(t);

                // Cell corners in physical space
                const float x0 = static_cast<float>(i)     * dx;
                const float x1 = static_cast<float>(i + 1) * dx;
                const float y0 = static_cast<float>(j)     * dy;
                const float y1 = static_cast<float>(j + 1) * dy;

                const std::size_t base = static_cast<std::size_t>((i * Ny + j) * 6);
                // Triangle 1: BL, BR, TR
                m_heatmap[base + 0] = {{x0, y0}, col};
                m_heatmap[base + 1] = {{x1, y0}, col};
                m_heatmap[base + 2] = {{x1, y1}, col};
                // Triangle 2: BL, TR, TL
                m_heatmap[base + 3] = {{x0, y0}, col};
                m_heatmap[base + 4] = {{x1, y1}, col};
                m_heatmap[base + 5] = {{x0, y1}, col};
            }
        }
    }

    void rebuildArrows()
    {
        const int Nx = m_snap.Nx, Ny = m_snap.Ny;
        const float dx = m_snap.Lx / static_cast<float>(Nx);
        const float dy = m_snap.Ly / static_cast<float>(Ny);

        // Draw arrows at every (stepI, stepJ) cell to avoid clutter
        const int stepI = std::max(1, Nx / 16);
        const int stepJ = std::max(1, Ny / 8);

        std::vector<sf::Vertex> verts;
        verts.reserve(static_cast<std::size_t>(Nx / stepI * Ny / stepJ * 4));

        // Compute max speed for scaling
        float velMax = 1e-6f;
        for (int k = 0; k < Nx * Ny; ++k) {
            const float s = std::sqrt(
                m_snap.velX[static_cast<std::size_t>(k)] * m_snap.velX[static_cast<std::size_t>(k)] +
                m_snap.velY[static_cast<std::size_t>(k)] * m_snap.velY[static_cast<std::size_t>(k)]);
            velMax = std::max(velMax, s);
        }

        const float arrowScale = 0.5f * std::min(dx * stepI, dy * stepJ) / velMax;

        for (int i = 0; i < Nx; i += stepI) {
            for (int j = 0; j < Ny; j += stepJ) {
                const std::size_t k = static_cast<std::size_t>(i * Ny + j);
                const float ux = m_snap.velX[k];
                const float uy = m_snap.velY[k];
                const float spd = std::sqrt(ux * ux + uy * uy);
                const float t = spd / velMax;
                const sf::Color col{255, 255, 255,
                    static_cast<uint8_t>(80 + static_cast<int>(120 * t))};

                // Centre of cell
                const float cx = (static_cast<float>(i) + 0.5f) * dx;
                const float cy = (static_cast<float>(j) + 0.5f) * dy;

                const float ex = cx + ux * arrowScale;
                const float ey = cy + uy * arrowScale;

                // Shaft
                verts.push_back({{cx, cy}, col});
                verts.push_back({{ex, ey}, col});

                // Arrowhead (two small lines from tip)
                const float len = spd * arrowScale;
                if (len > 1e-4f) {
                    const float nx = -uy / spd;
                    const float ny =  ux / spd;
                    const float hw = len * 0.25f;
                    verts.push_back({{ex, ey}, col});
                    verts.push_back({{ex - ux/spd*hw*1.2f + nx*hw*0.5f,
                                      ey - uy/spd*hw*1.2f + ny*hw*0.5f}, col});
                    verts.push_back({{ex, ey}, col});
                    verts.push_back({{ex - ux/spd*hw*1.2f - nx*hw*0.5f,
                                      ey - uy/spd*hw*1.2f - ny*hw*0.5f}, col});
                }
            }
        }

        m_arrows = sf::VertexArray(sf::PrimitiveType::Lines, verts.size());
        for (std::size_t i = 0; i < verts.size(); ++i)
            m_arrows[i] = verts[i];
    }

    void rebuildStreamlines()
    {
        m_streamlines.clear();

        const int Nx = m_snap.Nx, Ny = m_snap.Ny;

        // Max speed for coloring
        float velMax = 1e-6f;
        for (int k = 0; k < Nx * Ny; ++k) {
            const float s = std::sqrt(
                m_snap.velX[static_cast<std::size_t>(k)] * m_snap.velX[static_cast<std::size_t>(k)] +
                m_snap.velY[static_cast<std::size_t>(k)] * m_snap.velY[static_cast<std::size_t>(k)]);
            velMax = std::max(velMax, s);
        }

        const float dx = m_snap.Lx / static_cast<float>(Nx);
        const float dy = m_snap.Ly / static_cast<float>(Ny);

        // Seed along a grid: N_SX columns, N_SY rows
        constexpr int N_SX = 5;
        constexpr int N_SY = 8;

        for (int si = 0; si < N_SX; ++si) {
            for (int sj = 0; sj < N_SY; ++sj) {
                // Avoid seeding exactly on boundary
                const float x0 = dx * 0.6f + static_cast<float>(si)
                                * (m_snap.Lx - 1.2f * dx) / static_cast<float>(N_SX - 1);
                const float y0 = dy * 0.6f + static_cast<float>(sj)
                                * (m_snap.Ly - 1.2f * dy) / static_cast<float>(N_SY - 1);
                auto sl = integrateStreamline(x0, y0, velMax);
                if (sl.size() > 4)
                    m_streamlines.push_back(std::move(sl));
            }
        }
    }

    // Apply physical→screen transform and draw vertex array.
    // Called via the outer draw() to map the stored physical-coord geometry.
    // We override the standard draw() to inject a transform.
    //
    // The heatmap and arrows are stored in physical coords. In draw() we
    // compute the sf::Transform and pass it as sf::RenderStates.
    //
    // This method is called from draw() with the already-computed transform.
    friend class FieldViewWidget;

public:
    // Override: draw with an externally computed transform (physical→screen).
    // Called by the outer draw() method.
    void drawTransformed(sf::RenderTarget& target,
                         sf::FloatRect inner,
                         float Lx, float Ly) const
    {
        // Build affine transform: physical (x,y) → screen (sx,sy)
        // sx = inner.left + x * (inner.width  / Lx)
        // sy = inner.top  + inner.height - y * (inner.height / Ly)   [y-flip]
        //
        // sf::Transform is column-major 3x3:
        //   [a  b  tx]   a = sx/px,  b = 0
        //   [c  d  ty]   c = 0,      d = -sy/py  (negative for y-flip)
        //   [0  0   1]   tx = inner.left, ty = inner.top + inner.height
        const float a  = inner.size.x / Lx;
        const float d  = -inner.size.y / Ly;   // negative → y-flip
        const float tx = inner.position.x;
        const float ty = inner.position.y + inner.size.y;

        sf::Transform T;
        T = sf::Transform(a,  0.f, tx,
                          0.f, d,  ty,
                          0.f, 0.f, 1.f);

        sf::RenderStates rs;
        rs.transform = T;

        target.draw(m_heatmap, rs);

        if (m_mode != DisplayMode::Streamlines) {
            target.draw(m_arrows, rs);
        } else {
            for (const auto& sl : m_streamlines) {
                if (sl.size() < 2) continue;
                sf::VertexArray va(sf::PrimitiveType::LineStrip, sl.size());
                for (std::size_t i = 0; i < sl.size(); ++i)
                    va[i] = {{sl[i].x, sl[i].y}, sl[i].color};
                target.draw(va, rs);
            }
        }
    }
};

// ── FieldViewWidget ───────────────────────────────────────────────────────────
// Wraps FieldView and handles the full draw call including panel + colorbar.
// (Refactored out of FieldView::draw so draw() can use drawTransformed.)
class FieldViewWidget
{
public:
    explicit FieldViewWidget(const sf::Font& font) : m_view(font), m_font(font) {}

    void setMode(DisplayMode mode) { m_view.setMode(mode); }
    DisplayMode mode() const { return m_view.mode(); }
    void setSnapshot(const FieldSnapshot& snap) { m_view.setSnapshot(snap); }

    void draw(sf::RenderTarget& target, sf::FloatRect bounds) const
    {
        // Panel background
        RoundedPanel bg;
        bg.bounds       = bounds;
        bg.radius       = UILayout::RadiusPanel;
        bg.arcPts       = UILayout::ArcPts;
        bg.fill         = UIColors::Surface;
        bg.borderColor  = UIColors::Border;
        bg.borderThick  = UILayout::BorderThick;
        bg.shadowColor  = UIColors::Shadow;
        bg.shadowOffset = UILayout::ShadowOff;
        bg.draw(target);

        // Title
        const char* titles[] = {"Pressure", "Velocity Magnitude", "Streamlines"};
        sf::Text title(m_font, titles[static_cast<int>(m_view.mode())], 13u);
        title.setFillColor(UIColors::TextMuted);
        title.setPosition({bounds.position.x + UILayout::Pad,
                           bounds.position.y + UILayout::Pad});
        target.draw(title);

        const bool valid = m_view.m_snap.valid;

        if (!valid) {
            sf::Text hint(m_font, "Run a case to see the flow field", 13u);
            hint.setFillColor(UIColors::TextMuted);
            const auto hb = hint.getLocalBounds();
            hint.setPosition({
                bounds.position.x + (bounds.size.x - hb.size.x) * 0.5f,
                bounds.position.y + bounds.size.y * 0.5f});
            target.draw(hint);
            return;
        }

        constexpr float TitleH  = 24.f;
        constexpr float CbarW   = 28.f;
        constexpr float CbarGap = 8.f;
        const sf::FloatRect inner{
            {bounds.position.x + UILayout::Pad,
             bounds.position.y + UILayout::Pad + TitleH},
            {bounds.size.x - 2.f * UILayout::Pad - CbarW - CbarGap,
             bounds.size.y - 2.f * UILayout::Pad - TitleH}};

        const float Lx = m_view.m_snap.Lx;
        const float Ly = m_view.m_snap.Ly;

        // Clip drawing to inner area using a scissor-like approach:
        // SFML doesn't have built-in scissor, but RenderTexture would be overkill.
        // Just draw and rely on the panel background occluding anything outside.
        m_view.drawTransformed(target, inner, Lx, Ly);

        // Border around field area (screen-space)
        {
            const float fx = inner.position.x;
            const float fy = inner.position.y;
            const float fw = inner.size.x;
            const float fh = inner.size.y;
            sf::VertexArray border(sf::PrimitiveType::LineStrip, 5);
            border[0] = {{fx,      fy},      UIColors::Border};
            border[1] = {{fx + fw, fy},      UIColors::Border};
            border[2] = {{fx + fw, fy + fh}, UIColors::Border};
            border[3] = {{fx,      fy + fh}, UIColors::Border};
            border[4] = {{fx,      fy},      UIColors::Border};
            target.draw(border);
        }

        // Colorbar
        {
            const float cbX = inner.position.x + inner.size.x + CbarGap;
            const float cbY = inner.position.y;
            const float cbH = inner.size.y;
            constexpr int N = 64;
            sf::VertexArray cb(sf::PrimitiveType::TriangleStrip,
                               static_cast<std::size_t>(2 * N));
            for (int k = 0; k < N; ++k) {
                const float t  = static_cast<float>(k) / static_cast<float>(N - 1);
                const float cy = cbY + cbH * (1.f - t);
                const sf::Color col = (m_view.m_mode == DisplayMode::Pressure)
                                    ? Colormap::pressure(t)
                                    : Colormap::velocity(t);
                cb[static_cast<std::size_t>(2 * k)]     = {{cbX,           cy}, col};
                cb[static_cast<std::size_t>(2 * k + 1)] = {{cbX + CbarW,   cy}, col};
            }
            target.draw(cb);

            sf::Text tmax(m_font, fmtVal(m_view.m_scalarMax), 10u);
            tmax.setFillColor(UIColors::TextMuted);
            tmax.setPosition({cbX, cbY - 1.f});
            target.draw(tmax);

            sf::Text tmin(m_font, fmtVal(m_view.m_scalarMin), 10u);
            tmin.setFillColor(UIColors::TextMuted);
            const auto tmb = tmin.getLocalBounds();
            tmin.setPosition({cbX, cbY + cbH - tmb.size.y - tmb.position.y});
            target.draw(tmin);
        }
    }

private:
    static std::string fmtVal(float v)
    {
        char buf[32];
        if (std::abs(v) < 1e-3f || std::abs(v) > 999.f)
            std::snprintf(buf, sizeof(buf), "%.1e", static_cast<double>(v));
        else
            std::snprintf(buf, sizeof(buf), "%.3f", static_cast<double>(v));
        return buf;
    }

    FieldView      m_view;
    const sf::Font& m_font;
};
