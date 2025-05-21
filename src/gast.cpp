// $Source$
//--------------------------------------------------
// GAST
//--------------------------------------------------
// Proyecto_TT1: Greenwich Apparent Sidereal Time Implementation
//
// Copyright (c) 2025, Meysam Mahooti
//
// Created: 2025/05/08
//
/** @file GAST.cpp
 *  @brief Greenwich Apparent Sidereal Time calculation
 *
 *  @author Noel Del Rio Gonzalez
 *  @bug No known bugs
 */

#include "..\include\GAST.hpp"
#include "..\include\gmst.hpp"
#include "..\include\EqnEquinox.hpp"
#include "..\include\SAT_Const.hpp"
#include <cmath>

using namespace std;

double GAST(double Mjd_UT1) {
    // Calculate GAST: GMST + Equation of the Equinoxes, modulo 2pi
    double gstime = fmod(GMST(Mjd_UT1) + EqnEquinox(Mjd_UT1), 2.0 * Const::pi);
    
    // Ensure result is in [0, 2pi)
    if (gstime < 0.0) {
        gstime += 2.0 * Const::pi;
    }
    
    return gstime;
}