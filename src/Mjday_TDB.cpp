// $Source$
//--------------------------------------------------
// Mjday_TDB
//--------------------------------------------------
// Proyecto_TT1: Modified Julian Date (TDB) Calculation Implementation
//
// Copyright (c) 2025, Meysam Mahooti
//
// Created: 2025/04/23
//
/** @file Mjday_TDB.cpp
 *  @brief Implementation of Modified Julian Date (TDB) calculation
 *
 *  @author Noel Del Rio Gonzalez
 *  @bug No known bugs
 */

#include "Mjday_TDB.hpp"
#include "SAT_Const.hpp"

using namespace std;

double Mjday_TDB(double Mjd_TT) {
    double T_TT = (Mjd_TT - Const::MJD_J2000) / 36525.0;
    return Mjd_TT + (0.001658 * sin(628.3076 * T_TT + 6.2401)
                   + 0.000022 * sin(575.3385 * T_TT + 4.2970)
                   + 0.000014 * sin(1256.6152 * T_TT + 6.1969)
                   + 0.000005 * sin(606.9777 * T_TT + 4.0212)
                   + 0.000005 * sin(52.9691 * T_TT + 0.4444)
                   + 0.000002 * sin(21.3299 * T_TT + 5.5431)
                   + 0.000010 * sin(628.3076 * T_TT + 4.2490)) / 86400.0;
}