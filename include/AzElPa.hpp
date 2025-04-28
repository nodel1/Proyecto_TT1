// $Header$
//--------------------------------------------------
// AzElPa
//--------------------------------------------------
// Proyecto_TT1: Azimuth, Elevation, and Partials Calculation Header
//
// Copyright (c) 2025, Meysam Mahooti
//
// Created: 2025/04/23
//
/** @file AzElPa.hpp
 *  @brief Computes azimuth, elevation, and partials from local tangent coordinates
 *
 *  @author Noel Del Rio Gonzalez
 *  @bug No known bugs
 */

#ifndef _AZELPA_
#define _AZELPA_

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include "matrix.hpp"
#include "SAT_Const.hpp"

using namespace std;

//--------------------------------------------------
// AzElPa (const Matrix& s, double& Az, double& El, Matrix& dAds, Matrix& dEds)
//--------------------------------------------------
/**
 * @brief Computes azimuth, elevation, and partials from local tangent coordinates
 *
 * @param[in] s Topocentric local tangent coordinates (3x1, East-North-Zenith frame)
 * @param[out] Az Azimuth [rad]
 * @param[out] El Elevation [rad]
 * @param[out] dAds Partials of azimuth w.r.t. s (1x3)
 * @param[out] dEds Partials of elevation w.r.t. s (1x3)
 */
void AzElPa(const Matrix& s, double& Az, double& El, Matrix& dAds, Matrix& dEds);

#endif