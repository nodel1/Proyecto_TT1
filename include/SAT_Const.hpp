// $Header$
//--------------------------------------------------
// SAT_Const
//--------------------------------------------------
// Proyecto_TT1: Astronomical and Mathematical Constants Header
//
// Copyright (c) 2025, Meysam Mahooti
//
// Created: 2025/04/23
//
/** @file sat_const.hpp
 *  @brief Definition of astronomical and mathematical constants
 *
 *  @author Noel Del Rio Gonzalez
 *  @bug No known bugs
 */

#ifndef _SAT_CONST_
#define _SAT_CONST_

#include <cmath>

class Const {
public:
    // Mathematical constants
    static const double pi2;       // 2pi
    static const double Rad;       // Radians per degree
    static const double Deg;       // Degrees per radian
    static const double Arcs;      // Arcseconds per radian

    // General
    static const double MJD_J2000; // Modified Julian Date of J2000
    static const double T_B1950;   // Epoch B1950
    static const double c_light;   // Speed of light [m/s]; DE430
    static const double AU;        // Astronomical unit [m]; DE430

    // Physical parameters of the Earth, Sun, and Moon
    static const double R_Earth;   // Earth's radius [m]; DE430
    static const double f_Earth;   // Flattening; WGS-84
    static const double R_Sun;     // Sun's radius [m]; DE430
    static const double R_Moon;    // Moon's radius [m]; DE430

    // Earth rotation
    static const double omega_Earth; // [rad/s]; WGS-84

    // Gravitational coefficients
    static const double GM_Earth;    // [m^3/s^2]; DE430
    static const double GM_Sun;      // [m^3/s^2]; DE430
    static const double GM_Moon;     // [m^3/s^2]; DE430
    static const double GM_Mercury;  // [m^3/s^2]; DE430
    static const double GM_Venus;    // [m^3/s^2]; DE430
    static const double GM_Mars;     // [m^3/s^2]; DE430
    static const double GM_Jupiter;  // [m^3/s^2]; DE430
    static const double GM_Saturn;   // [m^3/s^2]; DE430
    static const double GM_Uranus;   // [m^3/s^2]; DE430
    static const double GM_Neptune;  // [m^3/s^2]; DE430
    static const double GM_Pluto;    // [m^3/s^2]; DE430

    // Solar radiation pressure
    static const double P_Sol;       // [N/m^2]; IERS 96
};

#endif