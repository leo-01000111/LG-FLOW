#pragma once
#include <SFML/Graphics/Color.hpp>
#include <algorithm>
#include <cstdint>

enum class DisplayMode { Pressure, Velocity, Streamlines };

namespace Colormap
{

inline sf::Color lerp(sf::Color a, sf::Color b, float t)
{
    return sf::Color{
        static_cast<uint8_t>(a.r + static_cast<int>((b.r - a.r) * t)),
        static_cast<uint8_t>(a.g + static_cast<int>((b.g - a.g) * t)),
        static_cast<uint8_t>(a.b + static_cast<int>((b.b - a.b) * t)),
        255u};
}

// Cool-warm diverging: blue → light-grey → red
inline sf::Color pressure(float t)
{
    t = std::clamp(t, 0.f, 1.f);
    if (t < 0.5f)
        return lerp(sf::Color{ 59,  76, 192}, sf::Color{220, 220, 220}, t * 2.f);
    return     lerp(sf::Color{220, 220, 220}, sf::Color{180,   4,  38}, (t - 0.5f) * 2.f);
}

// Hot: dark-blue → red → yellow → near-white
inline sf::Color velocity(float t)
{
    t = std::clamp(t, 0.f, 1.f);
    if (t < 0.33f)
        return lerp(sf::Color{ 10,  10,  60}, sf::Color{200,   0,   0}, t / 0.33f);
    if (t < 0.67f)
        return lerp(sf::Color{200,   0,   0}, sf::Color{255, 200,   0}, (t - 0.33f) / 0.34f);
    return     lerp(sf::Color{255, 200,   0}, sf::Color{255, 255, 220}, (t - 0.67f) / 0.33f);
}

} // namespace Colormap
