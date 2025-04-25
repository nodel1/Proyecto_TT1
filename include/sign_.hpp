// $Header$
//--------------------------------------------------
// sign_
//--------------------------------------------------
// Proyecto_TT1: Sign Function Header
//
// Copyright (c) 2025, Meysam Mahooti
//
// Created: 2025/04/23
//
/** @file sign_.hpp
 *  @brief Computes the absolute value of a with the sign of b
 *
 *  @author Noel Del Rio Gonzalez
 *  @bug No known bugs
 */

#ifndef _SIGN_
#define _SIGN_

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>

using namespace std;

//--------------------------------------------------
// sign_ (double a, double b)
//--------------------------------------------------
/**
 * @brief Computes the absolute value of a with the sign of b
 *
 * @param[in] a Value whose absolute value is taken
 * @param[in] b Value determining the sign of the result
 * @return double Absolute value of a with the sign of b
 */
double sign_(double a, double b);

#endif