// $Header$
//--------------------------------------------------
// R_x
//--------------------------------------------------
// Proyecto_TT1: Rotation Matrix around X-axis Header
//
// Copyright (c) 2025, Meysam Mahooti
//
// Created: 2025/04/23
//
/** @file R_x.hpp
 *  @brief Computes rotation matrix around the X-axis
 *
 *  @author Noel Del Rio Gonzalez
 *  @bug No known bugs
 */

#ifndef _R_X_
#define _R_X_

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include "matrix.hpp"

using namespace std;

//--------------------------------------------------
// R_x (double angle)
//--------------------------------------------------
/**
 * @brief Computes the 3x3 rotation matrix around the X-axis
 *
 * @param[in] angle Angle of rotation [rad]
 * @return Matrix 3x3 rotation matrix
 */
Matrix R_x(double angle);

#endif