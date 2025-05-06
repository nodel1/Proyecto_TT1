// $Header$
//--------------------------------------------------
// Position
//--------------------------------------------------
// Proyecto_TT1: Position Vector Calculation Header
//
// Copyright (c) 2025, Meysam Mahooti
//
// Created: 2025/04/23
//
/** @file Position.hpp
 *  @brief Computes position vector from geodetic coordinates
 *
 *  @author Noel Del Rio Gonzalez
 *  @bug No known bugs
 */

#ifndef _POSITION_
#define _POSITION_

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include "matrix.hpp"

using namespace std;

//--------------------------------------------------
// Position (double lon, double lat, double h, double r[3])
//--------------------------------------------------
/**
 * @brief Computes position vector from geodetic coordinates
 *
 * @param[in] lon Longitude [rad]
 * @param[in] lat Geodetic latitude [rad]
 * @param[in] h Altitude [m]
 */
 
Matrix Position(double lon,double  lat, double h);

#endif