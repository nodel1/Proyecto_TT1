// $Source$
//--------------------------------------------------
// LTC
//--------------------------------------------------
// Proyecto_TT1: Local Tangent Coordinates Transformation Implementation
//
// Copyright (c) 2025, Meysam Mahooti
//
// Created: 2025/05/08
//
/** @file LTC.cpp
 *  @brief Transformation from Greenwich meridian system to local tangent coordinates (East-North-Zenith)
 *
 *  @author Noel Del Rio Gonzalez
 *  @bug No known bugs
 */

#include "..\include\LTC.hpp"
#include "..\include\R_y.hpp"
#include "..\include\R_z.hpp"

using namespace std;

Matrix LTC(double lon, double lat) {
    // Compute rotation matrices
    Matrix Rz = R_z(lon);        // Store R_z(lon) in a variable
    Matrix Ry = R_y(-lat);       // Store R_y(-lat) in a variable
    Matrix M = Ry * Rz;          // Perform multiplication

    // Perform cyclic permutation of rows (1->2, 2->3, 3->1)
    for (int j = 1; j <= 3; j++) {
        double Aux = M(1, j);
        M(1, j) = M(2, j);
        M(2, j) = M(3, j);
        M(3, j) = Aux;
    }

    return M;
}