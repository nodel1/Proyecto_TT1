// $Header$
//--------------------------------------------------
// PrecMatrix
//--------------------------------------------------
// Proyecto_TT1: Precession Matrix Transformation Header
//
// Copyright (c) 2025, Meysam Mahooti
//
// Created: 2025/05/08
//
/** @file PrecMatrix.hpp
 *  @brief Header for precession transformation of equatorial coordinates
 *
 *  @author Noel Del Rio Gonzalez
 *  @bug No known bugs
 */

#ifndef PRECMATRIX_HPP
#define PRECMATRIX_HPP

#include "matrix.hpp"

Matrix PrecMatrix(double Mjd_1, double Mjd_2);

#endif