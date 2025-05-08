// $Source$
//--------------------------------------------------
// NutMatrix
//--------------------------------------------------
// Proyecto_TT1: Nutation Matrix Transformation Implementation
//
// Copyright (c) 2025, Meysam Mahooti
//
// Created: 2025/05/08
//
/** @file NutMatrix.cpp
 *  @brief Transformation from mean to true equator and equinox
 *
 *  @author Noel Del Rio Gonzalez
 *  @bug No known bugs
 */

#include "..\include\NutMatrix.hpp"
#include "..\include\MeanObliquity.hpp"
#include "..\include\NutAngles.hpp"
#include "..\include\R_x.hpp"
#include "..\include\R_z.hpp"

using namespace std;

Matrix NutMatrix(double Mjd_TT) {
    // Mean obliquity of the ecliptic
    double eps = MeanObliquity(Mjd_TT);

    // Nutation in longitude and obliquity
    double dpsi, deps;
    NutAngles(Mjd_TT, dpsi, deps);

    // Transformation from mean to true equator and equinox
    Matrix Rx1 = R_x(-eps - deps); // R_x(-eps - deps)
    Matrix Rz = R_z(-dpsi);        // R_z(-dpsi)
    Matrix Rx2 = R_x(eps);         // R_x(+eps)
    Matrix NutMat = Rx1 * Rz;      // R_x(-eps - deps) * R_z(-dpsi)
    NutMat = NutMat * Rx2;         // R_x(-eps - deps) * R_z(-dpsi) * R_x(+eps)

    return NutMat;
}