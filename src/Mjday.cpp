// $Source$
//--------------------------------------------------
// Mjday
//--------------------------------------------------
// Proyecto_TT1: Modified Julian Date Calculation Implementation
//
// Copyright (c) 2025, Meysam Mahooti
//
// Created: 2025/04/23
//
/** @file Mjday.cpp
 *  @brief Implementation of Modified Julian Date calculation
 *
 *  @author Noel Del Rio Gonzalez
 *  @bug No known bugs
 */

#include "..\include\Mjday.hpp"

using namespace std;

double Mjday(int yr, int mon, int day, int hr, int min, double sec) {
    double jd = 367.0 * yr
              - floor((7 * (yr + floor((mon + 9) / 12.0))) * 0.25)
              + floor(275 * mon / 9.0)
              + day + 1721013.5
              + ((sec / 60.0 + min) / 60.0 + hr) / 24.0;
    return jd - 2400000.5;
}