// $Source$
//--------------------------------------------------
// SAT_Const
//--------------------------------------------------
// Proyecto_TT1: Astronomical and Mathematical Constants Implementation
//
// Copyright (c) 2025, Meysam Mahooti
//
// Created: 2025/04/23
//
/** @file sat_const.cpp
 *  @brief Implementation of astronomical and mathematical constants
 *
 *  @author Noel360
 *  @bug No known bugs
 */

#include "SAT_Const.hpp"

// Mathematical constants
const double Const::pi2 = 2.0 * M_PI;
const double Const::Rad = M_PI / 180.0;
const double Const::Deg = 180.0 / M_PI;
const double Const::Arcs = 3600.0 * 180.0 / M_PI;

// General
const double Const::MJD_J2000 = 51544.5;
const double Const::T_B1950 = -0.500002108;
const double Const::c_light = 299792458.0;
const double Const::AU = 149597870700.0;

// Physical parameters of the Earth, Sun, and Moon
const double Const::R_Earth = 6378.1363e3;
const double Const::f_Earth = 1.0 / 298.257223563;
const double Const::R_Sun = 696000e3;
const double Const::R_Moon = 1738e3;

// Earth rotation
const double Const::omega_Earth = (15.04106717866910 / 3600.0) * (M_PI / 180.0);

// Gravitational coefficients
const double Const::GM_Earth = 398600.435436e9;
const double Const::GM_Sun = 132712440041.939400e9;
const double Const::GM_Moon = Const::GM_Earth / 81.30056907419062;
const double Const::GM_Mercury = 22031.780000e9;
const double Const::GM_Venus = 324858.592000e9;
const double Const::GM_Mars = 42828.375214e9;
const double Const::GM_Jupiter = 126712764.800000e9;
const double Const::GM_Saturn = 37940585.200000e9;
const double Const::GM_Uranus = 5794548.600000e9;
const double Const::GM_Neptune = 6836527.100580e9;
const double Const::GM_Pluto = 977.0000000000009e9;

// Solar radiation pressure
const double Const::P_Sol = 1367.0 / Const::c_light;