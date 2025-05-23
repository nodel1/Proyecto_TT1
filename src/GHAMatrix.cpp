// $Source$
//--------------------------------------------------
// GHAMatrix
//--------------------------------------------------
// Proyecto_TT1: Transformation Matrix Implementation
//
// Copyright (c) 2025, Meysam Mahooti
//
// Created: 2025/05/23
//
/** @file GHAMatrix.cpp
 *  @brief Computes the transformation matrix from true equator and equinox to Earth equator and Greenwich meridian system
 *
 *  @author Noel Del Rio
 *  @bug No known bugs
 */

#include "../include/GHAMatrix.hpp"
#include "../include/R_z.hpp"
#include "../include/GAST.hpp"

Matrix GHAMatrix(double Mjd_UT1) {
    return R_z(GAST(Mjd_UT1));
}