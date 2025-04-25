// $Header$
//--------------------------------------------------
// R_y
//--------------------------------------------------
// Proyecto_TT1: Rotation Matrix around Y-axis Header
//
// Copyright (c) 2025, Meysam Mahooti
//
// Created: 2025/04/23
//
/** @file R_y.hpp
 *  @brief Computes rotation matrix around the Y-axis
 *
 *  @author Noel Del Rio Gonzalez
 *  @bug No known bugs
 */

#ifndef _R_Y_
#define _R_Y_

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include "matrix.hpp"

using namespace std;

//--------------------------------------------------
// R_y (double angle)
//--------------------------------------------------
/**
 * @brief Computes the 3x3 rotation matrix around the Y-axis
 *
 * @param[in] angle Angle of rotation [rad]
 * @return Matrix 3x3 rotation matrix
 */
Matrix R_y(double angle);

#endif