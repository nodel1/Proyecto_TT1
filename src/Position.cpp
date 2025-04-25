// $Source$
//--------------------------------------------------
// Position
//--------------------------------------------------
// Proyecto_TT1: Position Vector Calculation Implementation
//
// Copyright (c) 2025, Meysam Mahooti
//
// Created: 2025/04/23
//
/** @file Position.cpp
 *  @brief Implementation of position vector calculation from geodetic coordinates
 *
 *  @author Noel Del Rio Gonzalez
 *  @bug No known bugs
 */

#include "Position.hpp"
#include "SAT_Const.hpp"

using namespace std;

void Position(double lon, double lat, double h, double r[3]) {
    double R_equ = Const::R_Earth;
    double f = Const::f_Earth;

    double e2 = f * (2.0 - f); // Square of eccentricity
    double CosLat = cos(lat);
    double SinLat = sin(lat);

    // Position vector
    double N = R_equ / sqrt(1.0 - e2 * SinLat * SinLat);

    r[0] = (N + h) * CosLat * cos(lon);
    r[1] = (N + h) * CosLat * sin(lon);
    r[2] = ((1.0 - e2) * N + h) * SinLat;
}