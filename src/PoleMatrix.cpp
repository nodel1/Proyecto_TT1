// $Source$
//--------------------------------------------------
// PoleMatrix
//--------------------------------------------------
// Proyecto_TT1: Pole Matrix Transformation Implementation
//
// Copyright (c) 2025, Meysam Mahooti
//
// Created: 2025/05/08
//
/** @file PoleMatrix.cpp
 *  @brief Transformation from pseudo Earth-fixed to Earth-fixed coordinates
 *
 *  @author Noel Del Rio Gonzalez
 *  @bug No known bugs
 */

#include "..\include\PoleMatrix.hpp"
#include "..\include\R_y.hpp"
#include "..\include\R_x.hpp"

using namespace std;

Matrix PoleMatrix(double xp, double yp) {
    // Compute rotation matrices
    Matrix Ry = R_y(-xp);  // R_y(-xp)
    Matrix Rx = R_x(-yp);  // R_x(-yp)

    // Transformation from pseudo Earth-fixed to Earth-fixed coordinates
    Matrix PoleMat = Ry * Rx;  // R_y(-xp) * R_x(-yp)

    return PoleMat;
}