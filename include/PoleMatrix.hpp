// $Header$
//--------------------------------------------------
// PoleMatrix
//--------------------------------------------------
// Proyecto_TT1: Pole Matrix Transformation Header
//
// Copyright (c) 2025, Meysam Mahooti
//
// Created: 2025/05/08
//
/** @file PoleMatrix.hpp
 *  @brief Header for transformation from pseudo Earth-fixed to Earth-fixed coordinates
 *
 *  @author Noel Del Rio Gonzalez
 *  @bug No known bugs
 */

#ifndef POLEMATRIX_HPP
#define POLEMATRIX_HPP

#include "matrix.hpp"

Matrix PoleMatrix(double xp, double yp);

#endif