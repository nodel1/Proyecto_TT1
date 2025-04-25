// $Header$
//--------------------------------------------------
// timediff
//--------------------------------------------------
// Proyecto_TT1: Time Differences Calculation Header
//
// Copyright (c) 2025, Meysam Mahooti
//
// Created: 2025/04/23
//
/** @file timediff.hpp
 *  @brief Computes time differences between various time scales
 *
 *  @author Noel Del Rio Gonzalez
 *  @bug No known bugs
 */

#ifndef _TIMEDIFF_
#define _TIMEDIFF_

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>

using namespace std;

//--------------------------------------------------
// timediff (double UT1_UTC, double TAI_UTC, double& UT1_TAI, double& UTC_GPS, double& UT1_GPS, double& TT_UTC, double& GPS_UTC)
//--------------------------------------------------
/**
 * @brief Computes time differences between various time scales
 *
 * @param[in] UT1_UTC Time difference between UT1 and UTC [s]
 * @param[in] TAI_UTC Time difference between TAI and UTC [s]
 * @param[out] UT1_TAI Time difference between UT1 and TAI [s]
 * @param[out] UTC_GPS Time difference between UTC and GPS [s]
 * @param[out] UT1_GPS Time difference between UT1 and GPS [s]
 * @param[out] TT_UTC Time difference between TT and UTC [s]
 * @param[out] GPS_UTC Time difference between GPS and UTC [s]
 */
void timediff(double UT1_UTC, double TAI_UTC, double& UT1_TAI, double& UTC_GPS, double& UT1_GPS, double& TT_UTC, double& GPS_UTC);

#endif