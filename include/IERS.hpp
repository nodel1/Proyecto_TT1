// $Header$
//--------------------------------------------------
// IERS
//--------------------------------------------------
// Proyecto_TT1: IERS Time and Polar Motion Data Management Header
//
// Copyright (c) 2025, Meysam Mahooti
//
// Created: 2025/04/23
//
/** @file IERS.hpp
 *  @brief Manages IERS time and polar motion data
 *
 *  @author Noel Del Rio Gonzalez
 *  @bug No known bugs
 */

#ifndef _IERS_
#define _IERS_

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include "matrix.hpp"
#include "SAT_Const.hpp"

using namespace std;

//--------------------------------------------------
// IERS (const Matrix& eop, double Mjd_UTC, char interp, double& x_pole, ...)
//--------------------------------------------------
/**
 * @brief Manages IERS time and polar motion data with specified interpolation
 *
 * @param[in] eop Earth orientation parameters matrix (at least 13 rows)
 * @param[in] Mjd_UTC Modified Julian Date in UTC
 * @param[in] interp Interpolation method ('l' for linear, 'n' for none)
 * @param[out] x_pole Polar coordinate x [rad]
 * @param[out] y_pole Polar coordinate y [rad]
 * @param[out] UT1_UTC UT1-UTC time difference [s]
 * @param[out] LOD Length of day [s]
 * @param[out] dpsi Nutation correction in longitude [rad]
 * @param[out] deps Nutation correction in obliquity [rad]
 * @param[out] dx_pole Polar coordinate correction x [rad]
 * @param[out] dy_pole Polar coordinate correction y [rad]
 * @param[out] TAI_UTC TAI-UTC time difference [s]
 */
void IERS(const Matrix& eop, double Mjd_UTC, char interp, double& x_pole, double& y_pole, double& UT1_UTC, double& LOD, double& dpsi, double& deps, double& dx_pole, double& dy_pole, double& TAI_UTC);

//--------------------------------------------------
// IERS (const Matrix& eop, double Mjd_UTC, double& x_pole, ...)
//--------------------------------------------------
/**
 * @brief Manages IERS time and polar motion data without interpolation
 *
 * @param[in] eop Earth orientation parameters matrix (at least 13 rows)
 * @param[in] Mjd_UTC Modified Julian Date in UTC
 * @param[out] x_pole Polar coordinate x [rad]
 * @param[out] y_pole Polar coordinate y [rad]
 * @param[out] UT1_UTC UT1-UTC time difference [s]
 * @param[out] LOD Length of day [s]
 * @param[out] dpsi Nutation correction in longitude [rad]
 * @param[out] deps Nutation correction in obliquity [rad]
 * @param[out] dx_pole Polar coordinate correction x [rad]
 * @param[out] dy_pole Polar coordinate correction y [rad]
 * @param[out] TAI_UTC TAI-UTC time difference [s]
 */
void IERS(const Matrix& eop, double Mjd_UTC, double& x_pole, double& y_pole, double& UT1_UTC, double& LOD, double& dpsi, double& deps, double& dx_pole, double& dy_pole, double& TAI_UTC);

#endif