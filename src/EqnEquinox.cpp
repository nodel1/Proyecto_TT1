// $Source$
//--------------------------------------------------
// EqnEquinox
//--------------------------------------------------
// Proyecto_TT1: Equation of the Equinoxes Implementation
//
// Copyright (c) 2025, Meysam Mahooti
//
// Created: 2025/05/07
//
/** @file EqnEquinox.cpp
 *  @brief Computes the equation of the equinoxes
 *
 *  @author Noel Del Rio Gonzalez
 *  @bug No known bugs
 */

#include "..\include\EqnEquinox.hpp"
#include "..\include\SAT_Const.hpp"
#include <cmath>

double EqnEquinox(double Mjd_TT) {
    // Nutation in longitude and obliquity
    double dpsi, deps;
    NutAngles(Mjd_TT, dpsi, deps);

    // Equation of the equinoxes
    double EqE = dpsi * std::cos(MeanObliquity(Mjd_TT));

    return EqE;
}