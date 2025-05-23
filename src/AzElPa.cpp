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
    const double pi2 = Const::pi2;

    // Calculate rho (horizontal distance in East-North plane)
    double rho = std::sqrt(s(1,1) * s(1,1) + s(2,1) * s(2,1));

    // Angles
    Az = std::atan2(s(1,1), s(2,1));
    if (Az < 0.0) {
        Az += pi2;
    }
    El = std::atan(s(3,1) / rho);

    // Partials
    dAds = zeros(1, 3); // Initialize 1x3 matrix with zeros
    dAds(1,1) = s(2,1) / (rho * rho);
    dAds(1,2) = -s(1,1) / (rho * rho);
    dAds(1,3) = 0.0;

    dEds = zeros(1, 3); // Initialize 1x3 matrix with zeros
    Matrix s_non_const = s; // Create non-const copy for dot
    double denom = dot(s_non_const, s_non_const);
    dEds(1,1) = -s(1,1) * s(3,1) / rho / denom;
    dEds(1,2) = -s(2,1) * s(3,1) / rho / denom;
    dEds(1,3) = rho / denom;
}