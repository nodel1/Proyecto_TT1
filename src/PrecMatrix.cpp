// $Source$
//--------------------------------------------------
// PrecMatrix
//--------------------------------------------------
// Proyecto_TT1: Precession Matrix Transformation Implementation
//
// Copyright (c) 2025, Meysam Mahooti
//
// Created: 2025/05/08
//
/** @file PrecMatrix.cpp
 *  @brief Precession transformation of equatorial coordinates
 *
 *  @author Noel Del Rio Gonzalez
 *  @bug No known bugs
 */

#include "..\include\PrecMatrix.hpp"
#include "..\include\R_y.hpp"
#include "..\include\R_z.hpp"
#include "..\include\SAT_Const.hpp"

using namespace std;

Matrix PrecMatrix(double Mjd_1, double Mjd_2) {
    // Time intervals in Julian centuries
    double T = (Mjd_1 - Const::MJD_J2000) / 36525.0;
    double dT = (Mjd_2 - Mjd_1) / 36525.0;

    // Precession angles
    double zeta = ((2306.2181 + (1.39656 - 0.000139 * T) * T) +
                   ((0.30188 - 0.000344 * T) + 0.017998 * dT) * dT) * dT / Const::Arcs;
    double z = zeta + ((0.79280 + 0.000411 * T) + 0.000205 * dT) * dT * dT / Const::Arcs;
    double theta = ((2004.3109 - (0.85330 + 0.000217 * T) * T) -
                    ((0.42665 + 0.000217 * T) + 0.041833 * dT) * dT) * dT / Const::Arcs;

    // Precession matrix
    Matrix Rz1 = R_z(-z);        // R_z(-z)
    Matrix Ry = R_y(theta);      // R_y(theta)
    Matrix Rz2 = R_z(-zeta);     // R_z(-zeta)
    Matrix PrecMat = Rz1 * Ry;   // R_z(-z) * R_y(theta)
    PrecMat = PrecMat * Rz2;     // R_z(-z) * R_y(theta) * R_z(-zeta)

    return PrecMat;
}