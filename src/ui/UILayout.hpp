#pragma once

namespace UILayout {
    // Window
    constexpr float WindowW       = 1280.f;
    constexpr float WindowH       =  720.f;

    // Spacing
    constexpr float Pad           =  12.f;
    constexpr float PadLg         =  20.f;
    constexpr float Gap           =   8.f;

    // Corner radii
    constexpr float RadiusPanel   =  10.f;
    constexpr float RadiusCard    =   8.f;
    constexpr float RadiusBtn     =   6.f;
    constexpr float RadiusInner   =   4.f;   // RadiusCard - Gap/2

    // Panels (x, y, w, h)
    constexpr float SideW         = 240.f;
    constexpr float SideX         =  Pad;
    constexpr float SideY         =  Pad;
    constexpr float SideH         = WindowH - 2.f * Pad;

    constexpr float StatusH       =  56.f;
    constexpr float StatusX       = SideX + SideW + Gap;
    constexpr float StatusY       =  Pad;
    constexpr float StatusW       = WindowW - StatusX - Pad;

    constexpr float FieldX        = StatusX;
    constexpr float FieldY        = StatusY + StatusH + Gap;
    constexpr float FieldW        = StatusW;
    constexpr float FieldH        = 400.f;

    constexpr float PlotX         = StatusX;
    constexpr float PlotY         = FieldY + FieldH + Gap;
    constexpr float PlotW         = StatusW;
    constexpr float PlotH         = WindowH - PlotY - Pad;

    // Drop shadow offset
    constexpr float ShadowOff     =   5.f;
    constexpr float ShadowOffHov  =   9.f;

    // Border thickness
    constexpr float BorderThick   =   1.5f;

    // Button height inside side panel
    constexpr float BtnH          =  40.f;

    // Points per corner arc (quality knob)
    constexpr int   ArcPts        =  10;
}
