// $Header$
//--------------------------------------------------
// MeasUpdate
//--------------------------------------------------
// Proyecto_TT1: Kalman Filter Measurement Update Header
//
// Copyright (c) 2025, Meysam Mahooti
//
// Created: 2025/05/21
//
/** @file MeasUpdate.hpp
 *  @brief Kalman filter measurement update implementation
 *
 *  @author Noel Del Rio Gonzalez
 *  @bug No known bugs
 */

#ifndef _MEASUPDATE_
#define _MEASUPDATE_

#include "matrix.hpp"
#include "SAT_Const.hpp"

using namespace std;

//--------------------------------------------------
// MeasUpdate (Matrix& x, double z, double g, double s, Matrix& G, Matrix& P, int n, Matrix& K)
//--------------------------------------------------
/**
 * @brief Performs Kalman filter measurement update step for scalar measurements
 *
 * @param[in,out] x State vector (n x 1)
 * @param[in] z Measurement value (scalar)
 * @param[in] g Measurement model output (scalar)
 * @param[in] s Measurement noise standard deviation (scalar)
 * @param[in] G Measurement sensitivity matrix (1 x n)
 * @param[in,out] P State covariance matrix (n x n)
 * @param[in] n State dimension
 * @param[out] K Kalman gain matrix (n x 1)
 */
void MeasUpdate(Matrix& x, double z, double g, double s, 
                const Matrix& G, Matrix& P, int n, Matrix& K);

#endif