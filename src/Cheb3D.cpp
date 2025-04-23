// $Source$
//--------------------------------------------------
// Cheb3D
//--------------------------------------------------
// Proyecto_TT1: Chebyshev Approximation Implementation
//
// Copyright (c) 2025, Meysam Mahooti
//
// Created: 2025/04/23
//
/** @file cheb3d.cpp
 *  @brief 3D Chebyshev approximation implementation
 *
 *  @author Noel Del Rio Gonzalez
 *  @bug No known bugs
 */

#include "cheb3d.hpp"
#include <stdexcept>

Matrix& Cheb3D(double t, int N, double Ta, double Tb, Matrix& Cx, Matrix& Cy, Matrix& Cz) {
    // Check validity
    if ((t < Ta) || (Tb < t)) {
        throw std::runtime_error("ERROR: Time out of range in Cheb3D::Value");
    }

    // Clenshaw algorithm
    double tau = (2*t - Ta - Tb)/(Tb - Ta);  

    Matrix& f1 = zeros(1,3);
    Matrix& f2 = zeros(1,3);

    for (int i = N; i >= 2; i--) {
        Matrix& old_f1 = f1;
        f1 = (f1 * (2*tau)) - f2 + Matrix{{Cx(i-1,0), Cy(i-1,0), Cz(i-1,0)}};
        f2 = old_f1;
    }

    return (f1 * tau) - f2 + Matrix{{Cx(0,0), Cy(0,0), Cz(0,0)}};
}