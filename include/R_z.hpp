// $Header$
//--------------------------------------------------
// R_z
//--------------------------------------------------
// Proyecto_TT1: Rotation Matrix around Z-axis Header
//
// Copyright (c) 2025, Meysam Mahooti
//
// Created: 2025/04/23
//
/** @file R_z.hpp
 *  @brief Computes rotation matrix around the Z-axis
 *
 *  @author Noel Del Rio Gonzalez
 *  @bug No known bugs
 */

#ifndef _R_Z_
#define _R_Z_

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include "matrix.hpp"

using namespace std;

//--------------------------------------------------
// R_z (double angle)
//--------------------------------------------------
/**
 * @brief Computes the 3x3 rotation matrix around the Z-axis
 *
 * @param[in] angle Angle of rotation [rad]
 * @return Matrix 3x3 rotation matrix
 */
Matrix R_z(double angle);

#endif