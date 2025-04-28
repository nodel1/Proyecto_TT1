// $Source$
//--------------------------------------------------
// MeanObliquity
//--------------------------------------------------
// Proyecto_TT1: Mean Obliquity Calculation Implementation
//
// Copyright (c) 2025, Meysam Mahooti
//
// Created: 2025/04/23
//
/** @file mean_obliquity.cpp
 *  @brief Mean obliquity of the ecliptic implementation
 *
 *  @author Noel Del Rio Gonzalez
 *  @bug No known bugs
 */

#include "..\include\MeanObliquity.hpp"
#include "SAT_Const.hpp"

double MeanObliquity(double Mjd_TT) {
    // Centuries since J2000
    double T = (Mjd_TT - Const::MJD_J2000)/36525.0;
    
    // Mean obliquity [rad]
    return Const::Rad * (84381.448/3600.0 - (46.8150 + (0.00059 - 0.001813*T)*T)*T/3600.0);
}