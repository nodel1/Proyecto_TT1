// $Source$
//--------------------------------------------------
// GMST
//--------------------------------------------------
// Proyecto_TT1: Greenwich Mean Sidereal Time Implementation
//
// Copyright (c) 2025, Meysam Mahooti
//
// Created: 2025/05/08
//
/** @file GMST.cpp
 *  @brief Greenwich Mean Sidereal Time calculation
 *
 *  @author Noel Del Rio Gonzalez
 *  @bug No known bugs
 */

#include "..\include\gmst.hpp"
#include "..\include\SAT_Const.hpp"
#include <cmath>

using namespace std;

double GMST(double Mjd_UT1) {
    const double Secs = 86400.0; // Seconds per day

    // Time calculations
    double Mjd_0 = floor(Mjd_UT1);
    double UT1 = Secs * (Mjd_UT1 - Mjd_0); // [s]
    double T_0 = (Mjd_0 - Const::MJD_J2000) / 36525.0;
    double T = (Mjd_UT1 - Const::MJD_J2000) / 36525.0;

    // GMST in seconds
    double gmst = 24110.54841 + 8640184.812866 * T_0 + 1.002737909350795 * UT1 +
                  (0.093104 - 6.2e-6 * T) * T * T; // [s]

    // Convert to radians [0, 2pi)
    double gmstime = 2.0 * Const::pi * fmod(gmst / Secs, 1.0); // [rad]

    return gmstime;
}