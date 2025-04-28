// $Source$
//--------------------------------------------------
// R_y
//--------------------------------------------------
// Proyecto_TT1: Rotation Matrix around Y-axis Implementation
//
// Copyright (c) 2025, Meysam Mahooti
//
// Created: 2025/04/23
//
/** @file R_y.cpp
 *  @brief Implementation of rotation matrix around the Y-axis
 *
 *  @author Noel Del Rio Gonzalez
 *  @bug No known bugs
 */

#include "..\include\R_y.hpp"

using namespace std;

Matrix R_y(double angle) {
    double C = cos(angle);
    double S = sin(angle);

    Matrix rotmat = zeros(3, 3); // Initialize 3x3 matrix with zeros

    rotmat(0, 0) =  C;  rotmat(0, 1) = 0.0;  rotmat(0, 2) = -S;
    rotmat(1, 0) = 0.0; rotmat(1, 1) = 1.0;  rotmat(1, 2) = 0.0;
    rotmat(2, 0) =  S;  rotmat(2, 1) = 0.0;  rotmat(2, 2) =  C;

    return rotmat;
}