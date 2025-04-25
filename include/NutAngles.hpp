// $Header$
//--------------------------------------------------
// NutAngles
//--------------------------------------------------
// Proyecto_TT1: Nutation Angles Calculation Header
//
// Copyright (c) 2025, Meysam Mahooti
//
// Created: 2025/04/23
//
/** @file NutAngles.hpp
 *  @brief Computes nutation in longitude and obliquity
 *
 *  @author Noel Del Rio Gonzalez
 *  @bug No known bugs
 */

#ifndef _NUTANGLES_
#define _NUTANGLES_

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>

using namespace std;

//--------------------------------------------------
// NutAngles (double Mjd_TT, double& dpsi, double& deps)
//--------------------------------------------------
/**
 * @brief Computes nutation in longitude and obliquity
 *
 * @param[in] Mjd_TT Modified Julian Date (Terrestrial Time)
 * @param[out] dpsi Nutation in longitude [rad]
 * @param[out] deps Nutation in obliquity [rad]
 */
void NutAngles(double Mjd_TT, double& dpsi, double& deps);

#endif