// $Header$
//--------------------------------------------------
// Mjday_TDB
//--------------------------------------------------
// Proyecto_TT1: Modified Julian Date (TDB) Calculation Header
//
// Copyright (c) 2025, Meysam Mahooti
//
// Created: 2025/04/23
//
/** @file Mjday_TDB.hpp
 *  @brief Computes the Modified Julian Date for Barycentric Dynamical Time
 *
 *  @author Noel Del Rio Gonzalez
 *  @bug No known bugs
 */

#ifndef _MJDAY_TDB_
#define _MJDAY_TDB_

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>

using namespace std;

//--------------------------------------------------
// Mjday_TDB (double Mjd_TT)
//--------------------------------------------------
/**
 * @brief Computes the Modified Julian Date for Barycentric Dynamical Time
 *
 * @param[in] Mjd_TT Modified Julian Date (Terrestrial Time)
 * @return double Modified Julian Date (TDB)
 */
double Mjday_TDB(double Mjd_TT);

#endif