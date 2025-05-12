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
#include "global.hpp"

using namespace std;

//--------------------------------------------------
// IERS (with interpolation)
//--------------------------------------------------
/**
 * @brief Computes IERS Earth rotation parameters with optional interpolation
 *
 * @param[in] Mjd_UTC Modified Julian Date in UTC [days]
 * @param[in] interp Interpolation method ('l' for linear, 'n' for none)
 * @param[out] x_pole Pole coordinate [rad]
 * @param[out] y_pole Pole coordinate [rad]
 * @param[out] UT1_UTC UT1-UTC time difference [s]
 * @param[out] LOD Length of day [s]
 * @param[out] dpsi Nutation correction [rad]
 * @param[out] deps Nutation correction [rad]
 * @param[out] dx_pole Pole coordinate correction [rad]
 * @param[out] dy_pole Pole coordinate correction [rad]
 * @param[out] TAI_UTC TAI-UTC time difference [s]
 */
void IERS(double Mjd_UTC, char interp, double& x_pole, double& y_pole, double& UT1_UTC, double& LOD, double& dpsi, double& deps, double& dx_pole, double& dy_pole, double& TAI_UTC);

//--------------------------------------------------
// IERS (no interpolation)
//--------------------------------------------------
/**
 * @brief Computes IERS Earth rotation parameters without interpolation
 *
 * @param[in] Mjd_UTC Modified Julian Date in UTC [days]
 * @param[out] x_pole Pole coordinate [rad]
 * @param[out] y_pole Pole coordinate [rad]
 * @param[out] UT1_UTC UT1-UTC time difference [s]
 * @param[out] LOD Length of day [s]
 * @param[out] dpsi Nutation correction [rad]
 * @param[out] deps Nutation correction [rad]
 * @param[out] dx_pole Pole coordinate correction [rad]
 * @param[out] dy_pole Pole coordinate correction [rad]
 * @param[out] TAI_UTC TAI-UTC time difference [s]
 */
void IERS(double Mjd_UTC, double& x_pole, double& y_pole, double& UT1_UTC, double& LOD, double& dpsi, double& deps, double& dx_pole, double& dy_pole, double& TAI_UTC);

#endif