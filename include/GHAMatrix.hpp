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
/** @file GHAMatrix.hpp
 *  @brief Computes the transformation matrix from true equator and equinox to Earth equator and Greenwich meridian system
 *
 *  @author Noel Del Rio
 *  @bug No known bugs
 */

#ifndef GHAMATRIX_HPP
#define GHAMATRIX_HPP

#include "Matrix.hpp"

/**
 * @brief Computes the Greenwich Hour Angle matrix
 * @param Mjd_UT1 Modified Julian Date UT1
 * @return GHAmat Greenwich Hour Angle matrix (3x3)
 */
Matrix GHAMatrix(double Mjd_UT1);

#endif // GHAMATRIX_HPP