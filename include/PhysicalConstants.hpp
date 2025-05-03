//==============================================================
// CGSConstants.hpp
// Fundamental physical constants in CGS (centimeter-gram-second) units
// All constants are defined as float for consistency across code base
//==============================================================

#pragma once

struct CGSConstants {
    static constexpr float G     = 6.67430e-8f;      // g^-1 cm^3 s^-2 (Gravitational constant)
    static constexpr float Msun  = 1.98841e33f;      // g (Solar mass)
    static constexpr float AU    = 1.49598e13f;      // cm (Astronomical Unit)
    static constexpr float year  = 3.15576e7f;       // s (Julian year)
    static constexpr float c     = 2.99792e10f;      // cm/s (Speed of light)
    static constexpr float hbar  = 1.05457e-27f;     // g·cm²/s (Reduced Planck constant)
    static constexpr float kB    = 1.38065e-16f;     // g·cm²/K·s² (Boltzmann constant)
};