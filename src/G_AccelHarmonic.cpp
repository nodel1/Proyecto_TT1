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
    Matrix da = zeros(3, 1);
    Matrix r_plus;  // Declarada fuera del bucle
    Matrix r_minus; // Declarada fuera del bucle
    Matrix a_plus;  // Declarada fuera del bucle
    Matrix a_minus; // Declarada fuera del bucle

    for (int i = 1; i <= 3; ++i) {
        dr = zeros(3, 1);
        dr(i, 1) = d;
        r_plus = r + dr * 0.5;
        r_minus = r - dr * 0.5;
        

        Matrix temp_plus = AccelHarmonic(r_plus, U, n_max, m_max); // temp_plus es un lvalue
        a_plus = temp_plus; // Usa operator=(Matrix&), que funciona con lvalues
        
        Matrix temp_minus = AccelHarmonic(r_minus, U, n_max, m_max); // temp_minus es un lvalue
        a_minus = temp_minus; // Usa operator=(Matrix&), que funciona con lvalues
        
        da = a_plus - a_minus;
        for (int j = 1; j <= 3; ++j) {
            G(j, i) = da(j, 1) / d;
        }
    }

    return G;
}