// $Source$
//--------------------------------------------------
// AzElPa
//--------------------------------------------------
// Proyecto_TT1: Azimuth, Elevation, and Partials Calculation Implementation
//
// Copyright (c) 2025, Meysam Mahooti
//
// Created: 2025/04/23
//
/** @file AzElPa.cpp
 *  @brief Implementation of azimuth, elevation, and partials calculation
 *
 *  @author Noel Del Rio Gonzalez
 *  @bug No known bugs
 */

#include "..\include\AzElPa.hpp"

using namespace std;

void AzElPa(const Matrix& s, double& Az, double& El, Matrix& dAds, Matrix& dEds) {
    const double pi2 = 2.0 * Const::pi;

    // Calculate rho (horizontal distance in East-North plane)
    double rho = sqrt(s(0,0) * s(0,0) + s(1,0) * s(1,0));

    // Angles
    Az = atan2(s(0,0), s(1,0));
    if (Az < 0.0) {
        Az += pi2;
    }
    El = atan(s(2,0) / rho);

    // Partials
    dAds = Matrix(1, 3); // Initialize 1x3 matrix
    dAds(0,0) = s(1,0) / (rho * rho);
    dAds(0,1) = -s(0,0) / (rho * rho);
    dAds(0,2) = 0.0;

    dEds = Matrix(1, 3); // Initialize 1x3 matrix
    double denom = dot(s, s);
    dEds(0,0) = -s(0,0) * s(2,0) / rho / denom;
    dEds(0,1) = -s(1,0) * s(2,0) / rho / denom;
    dEds(0,2) = rho / denom;
}