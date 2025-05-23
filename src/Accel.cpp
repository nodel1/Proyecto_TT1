// $Source$
//--------------------------------------------------
// Accel
//--------------------------------------------------
// Proyecto_TT1: Satellite Acceleration Computation
//
// Copyright (c) 2025, Meysam Mahooti
//
// Created: 2025/05/23
//
/** @file Accel.cpp
 *  @brief Computes the acceleration of an Earth-orbiting satellite
 *
 *  @author Noel Del Rio
 *  @bug No known bugs
 */

#include "../include/Accel.hpp"
#include "../include/IERS.hpp"
#include "../include/timediff.hpp"
#include "../include/PrecMatrix.hpp"
#include "../include/NutMatrix.hpp"
#include "../include/PoleMatrix.hpp"
#include "../include/GHAMatrix.hpp"
#include "../include/Mjday_TDB.hpp"
#include "../include/JPL_Eph_DE430.hpp"
#include "../include/AccelHarmonic.hpp"
#include "../include/AccelPointMass.hpp"
#include "../include/global.hpp"

Matrix Accel(double x, Matrix Y) {


    double x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC;

    double Mjd_UTC = AuxParam.Mjd_UTC + x / 86400.0;
    IERS(Mjd_UTC, 'l', x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC);



	double UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC;
	timediff(UT1_UTC, TAI_UTC, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);


    double Mjd_UT1 = AuxParam.Mjd_UTC + x / 86400.0 + UT1_UTC / 86400.0;
    double Mjd_TT = AuxParam.Mjd_UTC + x / 86400.0 + TT_UTC / 86400.0;


    Matrix P = PrecMatrix(Const::MJD_J2000, Mjd_TT);
    Matrix N = NutMatrix(Mjd_TT);
    Matrix T = N * P;


	Matrix GHA = GHAMatrix(Mjd_UT1);
	Matrix temp = GHA * T;
	Matrix E = PoleMatrix(x_pole, y_pole) * temp;

    // Get planetary positions
    double MJD_TDB = Mjday_TDB(Mjd_TT);
    Matrix r_Mercury, r_Venus, r_Earth, r_Mars, r_Jupiter, r_Saturn, r_Uranus, r_Neptune, r_Pluto, r_Moon, r_Sun;
    JPL_Eph_DE430(MJD_TDB, r_Mercury, r_Venus, r_Earth, r_Mars, r_Jupiter, r_Saturn, r_Uranus, r_Neptune, r_Pluto, r_Moon, r_Sun);

    // Extract position vector
    Matrix r(3, 1);
    for (int i = 1; i <= 3; ++i) {
        r(i, 1) = Y(i, 1);
    }

    // Compute acceleration due to harmonic gravity field
    Matrix a = AccelHarmonic(r, E, AuxParam.n, AuxParam.m);

    // Luni-solar perturbations
    if (AuxParam.sun) {
        a = a + AccelPointMass(r, r_Sun, Const::GM_Sun);
    }
    if (AuxParam.moon) {
        a = a + AccelPointMass(r, r_Moon, Const::GM_Moon);
    }

    // Planetary perturbations
    if (AuxParam.planets) {
        a = a + AccelPointMass(r, r_Mercury, Const::GM_Mercury);
        a = a + AccelPointMass(r, r_Venus, Const::GM_Venus);
        a = a + AccelPointMass(r, r_Mars, Const::GM_Mars);
        a = a + AccelPointMass(r, r_Jupiter, Const::GM_Jupiter);
        a = a + AccelPointMass(r, r_Saturn, Const::GM_Saturn);
        a = a + AccelPointMass(r, r_Uranus, Const::GM_Uranus);
        a = a + AccelPointMass(r, r_Neptune, Const::GM_Neptune);
        a = a + AccelPointMass(r, r_Pluto, Const::GM_Pluto);
    }


    Matrix dY = zeros(6, 1);
    for (int i = 1; i <= 3; ++i) {
        dY(i, 1) = Y(i + 3, 1); // Velocity components
        dY(i + 3, 1) = a(i, 1); // Acceleration components
    }

    return dY;
}