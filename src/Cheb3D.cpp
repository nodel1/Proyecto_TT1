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

#include "..\include\Cheb3d.hpp"
#include <stdexcept>

Matrix Cheb3D(double t, int N, double Ta, double Tb, Matrix& Cx, Matrix& Cy, Matrix& Cz) {
    // Check validity
    if (t < Ta || t > Tb) {
		cout << "ERROR: Time out of range in Cheb3D::Value\n";
        exit(EXIT_FAILURE);
    }

    // Clenshaw algorithm
    double tau = (2.0 * t - Ta - Tb) / (Tb - Ta);

    Matrix f1 = zeros(1, 3);
    Matrix f2 = zeros(1, 3);

    for (int i = N; i >= 2; i--) {
        Matrix old_f1 = f1;
        // Compute 2*tau*f1 - f2 + [Cx(i), Cy(i), Cz(i)]
        Matrix temp(1, 3);
        temp(1, 1) = Cx(i, 1);
        temp(1, 2) = Cy(i, 1);
        temp(1, 3) = Cz(i, 1);
        f1 = (f1 * (2.0 * tau)) - f2 + temp;
        f2 = old_f1;
    }

    // Compute tau*f1 - f2 + [Cx(1), Cy(1), Cz(1)]
    Matrix temp(1, 3);
    temp(1, 1) = Cx(1, 1);
    temp(1, 2) = Cy(1, 1);
    temp(1, 3) = Cz(1, 1);
    return (f1 * tau) - f2 + temp;
}
