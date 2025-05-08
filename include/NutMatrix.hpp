// $Header$
//--------------------------------------------------
// NutMatrix
//--------------------------------------------------
// Proyecto_TT1: Nutation Matrix Transformation Header
//
// Copyright (c) 2025, Meysam Mahooti
//
// Created: 2025/05/08
//
/** @file NutMatrix.hpp
 *  @brief Header for transformation from mean to true equator and equinox
 *
 *  @author Noel Del Rio Gonzalez
 *  @bug No known bugs
 */

#ifndef NUTMATRIX_HPP
#define NUTMATRIX_HPP

#include "matrix.hpp"

Matrix NutMatrix(double Mjd_TT);

#endif