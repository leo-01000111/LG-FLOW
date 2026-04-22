#pragma once
#include <SFML/Graphics/Color.hpp>

namespace UIColors {
    constexpr sf::Color Background  { 18,  18,  24, 255};
    constexpr sf::Color Surface     { 32,  32,  42, 255};
    constexpr sf::Color SurfaceHigh { 44,  44,  58, 255};
    constexpr sf::Color Accent      { 94, 161, 235, 255};
    constexpr sf::Color AccentDim   { 60, 110, 180, 255};
    constexpr sf::Color Border      { 60,  60,  75, 255};
    constexpr sf::Color Shadow      {  0,   0,   0, 120};
    constexpr sf::Color TextPrimary {220, 220, 230, 255};
    constexpr sf::Color TextMuted   {140, 140, 155, 180};
    constexpr sf::Color Success     { 80, 200, 120, 255};
    constexpr sf::Color Warning     {230, 170,  60, 255};
    constexpr sf::Color Danger      {220,  70,  70, 255};
    constexpr sf::Color PlotLine    { 94, 161, 235, 255};
    constexpr sf::Color PlotLine2   {200, 120,  80, 255};
    constexpr sf::Color Transparent {  0,   0,   0,   0};
}
