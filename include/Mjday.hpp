// $Header$
//--------------------------------------------------
// Mjday
//--------------------------------------------------
// Proyecto_TT1: Modified Julian Date Calculation Header
//
// Copyright (c) 2025, Meysam Mahooti
//
// Created: 2025/04/23
//
/** @file Mjday.hpp
 *  @brief Computes the Modified Julian Date from date and time
 *
 *  @author Noel Del Rio Gonzalez
 *  @bug No known bugs
 */

#ifndef _MJDAY_
#define _MJDAY_

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>

using namespace std;

//--------------------------------------------------
// Mjday (int yr, int mon, int day, int hr, int min, double sec)
//--------------------------------------------------
/**
 * @brief Computes the Modified Julian Date from date and time components
 *
 * @param[in] yr Year
 * @param[in] mon Month
 * @param[in] day Day
 * @param[in] hr Universal time hour (default: 0)
 * @param[in] min Universal time minute (default: 0)
 * @param[in] sec Universal time second (default: 0.0)
 * @return double Modified Julian Date
 */
double Mjday(int yr, int mon, int day, int hr = 0, int min = 0, double sec = 0.0);

#endif