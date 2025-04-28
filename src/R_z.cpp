// $Source$
//--------------------------------------------------
// R_z
//--------------------------------------------------
// Proyecto_TT1: Rotation Matrix around Z-axis Implementation
//
// Copyright (c) 2025, Meysam Mahooti
//
// Created: 2025/04/23
//
/** @file R_z.cpp
 *  @brief Implementation of rotation matrix around the Z-axis
 *
 *  @author Noel Del Rio Gonzalez
 *  @bug No known bugs
 */

#include "..\include\R_z.hpp"

using namespace std;

Matrix R_z(double angle) {
    double C = cos(angle);
    double S = sin(angle);

    Matrix rotmat = zeros(3, 3); // Initialize 3x3 matrix with zeros

    rotmat(0, 0) =  C;  rotmat(0, 1) =  S;  rotmat(0, 2) = 0.0;
    rotmat(1, 0) = -S;  rotmat(1, 1) =  C;  rotmat(1, 2) = 0.0;
    rotmat(2, 0) = 0.0; rotmat(2, 1) = 0.0; rotmat(2, 2) = 1.0;

    return rotmat;
}