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

    std::cout << "Entering Accel with x = " << x << std::endl;

    double x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC;
    double Mjd_UTC = AuxParam.Mjd_UTC + x / 86400.0;
    std::cout << "Calling IERS with Mjd_UTC = " << Mjd_UTC << std::endl;
    IERS(Mjd_UTC, 'l', x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC);

    std::cout << "despues de hacer iers" << std::endl;



    double UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC;
    std::cout << "Calling timediff..." << std::endl;
    timediff(UT1_UTC, TAI_UTC, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);

    double Mjd_UT1 = AuxParam.Mjd_UTC + x / 86400.0 + UT1_UTC / 86400.0;
    double Mjd_TT = AuxParam.Mjd_UTC + x / 86400.0 + TT_UTC / 86400.0;
    std::cout << "Mjd_UT1 = " << Mjd_UT1 << ", Mjd_TT = " << Mjd_TT << std::endl;

    std::cout << "Computing transformation matrices..." << std::endl;
    Matrix P = PrecMatrix(Const::MJD_J2000, Mjd_TT);
    Matrix N = NutMatrix(Mjd_TT);
    Matrix T = N * P;
    Matrix GHA = GHAMatrix(Mjd_UT1);
    Matrix temp = GHA * T;
    Matrix E = PoleMatrix(x_pole, y_pole) * temp;

    std::cout << "Calling JPL_Eph_DE430..." << std::endl;
    double MJD_TDB = Mjday_TDB(Mjd_TT);
    Matrix r_Mercury, r_Venus, r_Earth, r_Mars, r_Jupiter, r_Saturn, r_Uranus, r_Neptune, r_Pluto, r_Moon, r_Sun;
    JPL_Eph_DE430(MJD_TDB, r_Mercury, r_Venus, r_Earth, r_Mars, r_Jupiter, r_Saturn, r_Uranus, r_Neptune, r_Pluto, r_Moon, r_Sun);

    std::cout << "Extracting position vector..." << std::endl;
    Matrix r(3, 1);
    for (int i = 1; i <= 3; ++i) {
        r(i, 1) = Y(i, 1);
    }

    std::cout << "Calling AccelHarmonic with n = " << AuxParam.n << ", m = " << AuxParam.m << std::endl;
    Matrix a = AccelHarmonic(r, E, AuxParam.n, AuxParam.m);


    std::cout << "MATRIZ a" << std::endl;
		std::cout << a;
		
		
    std::cout << "MATRIZ r" << std::endl;
		std::cout << r;
	
	
	    std::cout << "Matriz r_sun" << std::endl;
			std::cout << r_Sun;


    std::cout << "Adding luni-solar perturbations..." << std::endl;
    if (AuxParam.sun) {
        a = a + AccelPointMass(r, transpose(r_Sun), Const::GM_Sun);
    }
    if (AuxParam.moon) {
        a = a + AccelPointMass(r, transpose(r_Moon), Const::GM_Moon);
    }

    if (AuxParam.planets) {
        std::cout << "Adding planetary perturbations..." << std::endl;
        a = a + AccelPointMass(r, transpose(r_Mercury), Const::GM_Mercury);
        a = a + AccelPointMass(r, transpose(r_Venus), Const::GM_Venus);
        a = a + AccelPointMass(r, transpose(r_Mars), Const::GM_Mars);
        a = a + AccelPointMass(r, transpose(r_Jupiter), Const::GM_Jupiter);
        a = a + AccelPointMass(r, transpose(r_Saturn), Const::GM_Saturn);
        a = a + AccelPointMass(r, transpose(r_Uranus), Const::GM_Uranus);
        a = a + AccelPointMass(r, transpose(r_Neptune), Const::GM_Neptune);
        a = a + AccelPointMass(r, transpose(r_Pluto), Const::GM_Pluto);
    }


    std::cout << "MATRIZ a" << std::endl;
		std::cout << a;
		
		    std::cout << "MATRIZ Y" << std::endl;
		std::cout << Y;
		
		
		
		
    std::cout << "Constructing dY..." << std::endl;
    Matrix dY = zeros(6, 1);
	dY(1,1) = Y(4,1);
	dY(2,1) = Y(5,1);
	dY(3,1) = Y(6,1);
	dY(4,1) = a(1,1);
	dY(5,1) = a(2,1);
	dY(6,1) = a(3,1);
	
	

    std::cout << "Returning dY from Accel..." << std::endl;
    return dY;

}