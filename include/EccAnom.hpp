// $Header$
//--------------------------------------------------
// EccAnom
//--------------------------------------------------
// Proyecto_TT1: Eccentric Anomaly Calculation Header
//
// Copyright (c) 2025, Meysam Mahooti
//
// Created: 2025/04/23
//
/** @file ecc_anom.hpp
 *  @brief Eccentric anomaly computation for elliptic orbits
 *
 *  @author Noel Del Rio Gonzalez
 *  @bug No known bugs
 */
 
#ifndef _ECCANOM_
#define _ECCANOM_
 
#include "..\include\SAT_Const.hpp"
#include <cmath>
#include <iostream>

using namespace std;

//--------------------------------------------------
// EccAnom (double M, double e)
//--------------------------------------------------
/**
 * @brief Computes the eccentric anomaly for elliptic orbits
 *
 * @param[in] M Mean anomaly in [rad]
 * @param[in] e Eccentricity of the orbit [0,1]
 * @return double Eccentric anomaly in [rad]
 *
 * @note Uses Newton-Raphson iteration method
 */
double EccAnom(double M, double e);

#endif