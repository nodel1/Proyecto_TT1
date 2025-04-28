// $Source$
//--------------------------------------------------
// Frac
//--------------------------------------------------
// Proyecto_TT1: Fractional Part Calculation Implementation
//
// Copyright (c) 2025, Meysam Mahooti
//
// Created: 2025/04/23
//
/** @file frac.cpp
 *  @brief Fractional part calculation implementation
 *
 *  @author Noel Del Rio Gonzalez
 *  @bug No known bugs
 */

#include "..\include\Frac.hpp"
#include <cmath>

double Frac(double x) {
    return x - floor(x);
}