// $Header$
//--------------------------------------------------
// Legendre
//--------------------------------------------------
// Proyecto_TT1: Associated Legendre Functions Header
//
// Copyright (c) 2025, Meysam Mahooti
//
// Created: 2025/04/23
//
/** @file Legendre.hpp
 *  @brief Computes normalized associated Legendre functions and their derivatives
 *
 *  @author Noel Del Rio Gonzalez
 *  @bug No known bugs
 */

#ifndef _LEGENDRE_
#define _LEGENDRE_

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <utility>
#include "matrix.hpp"

using namespace std;

//--------------------------------------------------
// Legendre (int n, int m, double fi)
//--------------------------------------------------
/**
 * @brief Computes normalized associated Legendre functions and their derivatives
 *
 * @param[in] n Maximum degree of Legendre functions
 * @param[in] m Maximum order of Legendre functions
 * @param[in] fi Angle [rad]
 */
void Legendre(int n, int m, double fi, Matrix& pnm, Matrix& dpnm);

#endif