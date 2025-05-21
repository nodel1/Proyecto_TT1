// $Source$
//--------------------------------------------------
// MeasUpdate
//--------------------------------------------------
// Proyecto_TT1: Kalman Filter Measurement Update Implementation
//
// Copyright (c) 2025, Meysam Mahooti
//
// Created: 2025/05/21
//
/** @file MeasUpdate.cpp
 *  @brief Implementation of Kalman filter measurement update step for scalar measurements
 *
 *  @author Noel Del Rio Gonzalez
 *  @bug No known bugs
 */

#include "..\include\MeasUpdate.hpp"
#include <cmath>

using namespace std;

void MeasUpdate(Matrix& x, double z, double g, double s, 
                const Matrix& G, Matrix& P, int n, Matrix& K)
{
    // Kalman gain calculation (scalar measurement version)
    Matrix G_transpose = G.transpose();
    double S = s*s + (G*P*G_transpose)(1,1);  // Innovation covariance (scalar)
    K = (P*G_transpose)*(1.0/S);  // Kalman gain (n x 1)

    // State update
    double innovation = z - g;
    x = x + K*innovation;

    // Covariance update (Joseph form for better numerical stability)
    Matrix I = eye(n);
    Matrix I_KG = I - K*G;
    P = I_KG*P*I_KG.transpose() + K*(s*s)*K.transpose();
}