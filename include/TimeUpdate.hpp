// $Header$
//--------------------------------------------------
// TimeUpdate
//--------------------------------------------------
// Proyecto_TT1: Covariance Matrix Update Header
//
// Copyright (c) 2025, Meysam Mahooti
//
// Created: 2025/04/23
//
/** @file TimeUpdate.hpp
 *  @brief Updates the covariance matrix using the state transition matrix and process noise
 *
 *  @author Noel Del Rio Gonzalez
 *  @bug No known bugs
 */

#ifndef _TIMEUPDATE_
#define _TIMEUPDATE_

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include "matrix.hpp"

using namespace std;

//--------------------------------------------------
// TimeUpdate (Matrix P, Matrix Phi, Matrix Qdt)
//--------------------------------------------------
/**
 * @brief Updates the covariance matrix using the state transition matrix and process noise
 *
 * @param[in] P Covariance matrix
 * @param[in] Phi State transition matrix
 * @param[in] Qdt Process noise matrix
 * @return Matrix Updated covariance matrix
 */
Matrix TimeUpdate(const Matrix& P, const Matrix& Phi, const Matrix& Qdt);

//--------------------------------------------------
// TimeUpdate (Matrix P, Matrix Phi)
//--------------------------------------------------
/**
 * @brief Updates the covariance matrix using the state transition matrix, with zero process noise
 *
 * @param[in] P Covariance matrix
 * @param[in] Phi State transition matrix
 * @return Matrix Updated covariance matrix
 */
Matrix TimeUpdate(const Matrix& P, const Matrix& Phi);

#endif