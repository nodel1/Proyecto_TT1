// $Source$
//--------------------------------------------------
// G_AccelHarmonic
//--------------------------------------------------
// Proyecto_TT1: Gradient of Harmonic Gravity Field Implementation
//
// Copyright (c) 2025, Meysam Mahooti
//
// Created: 2025/05/23
//
/** @file G_AccelHarmonic.cpp
 *  @brief Computes the gradient of the Earth's harmonic gravity field
 *
 *  @author Noel Del Rio Gonzalez
 *  @bug No known bugs
 */

#include "..\include\G_AccelHarmonic.hpp"
#include "..\include\AccelHarmonic.hpp"

Matrix G_AccelHarmonic(Matrix& r, Matrix& U, int n_max, int m_max) {
    const double d = 1.0;
    Matrix G = zeros(3, 3);
    Matrix dr = zeros(3, 1);

    for (int i = 1; i <= 3; ++i) {
        dr = zeros(3, 1);
        dr(i, 1) = d;
        Matrix r_plus = r + dr * 0.5;
        Matrix r_minus = r - dr * 0.5;
        Matrix a_plus = AccelHarmonic(r_plus, U, n_max, m_max);
        Matrix a_minus = AccelHarmonic(r_minus, U, n_max, m_max);
        Matrix da = a_plus - a_minus;
        for (int j = 1; j <= 3; ++j) {
            G(j, i) = da(j, 1) / d;
        }
    }

    return G;
}