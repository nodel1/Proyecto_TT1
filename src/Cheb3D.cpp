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
    // Comprobar validez de tiempo
    if (t < Ta || t > Tb) {
        std::cout << "ERROR: Time out of range in Cheb3D::Value\n";
        std::exit(EXIT_FAILURE);
    }

    // Comprobar tamaño de las matrices
    if (Cx.n_column != N || Cy.n_column != N || Cz.n_column != N) {
        throw std::runtime_error("Cheb3D: Incorrect size of Cx/Cy/Cz");
    }

    // Calcular tau
    double tau = (2.0 * t - Ta - Tb) / (Tb - Ta);

    // Inicializar f1 y f2
    Matrix f1 = zeros(1, 3);
    Matrix f2 = zeros(1, 3);
    Matrix old_f1 = zeros(1, 3);

    // Algoritmo de Clenshaw
    for (int i = N; i >= 2; i--) {
        old_f1 = f1;
        Matrix vaux = zeros(1, 3);
        vaux(1, 1) = Cx(1, i);
        vaux(1, 2) = Cy(1, i);
        vaux(1, 3) = Cz(1, i);
        f1 = (old_f1 * (2.0 * tau)) - f2 + vaux;
        f2 = old_f1;
    }

    // Calcular ChebApp
    Matrix vaux = zeros(1, 3);
    vaux(1, 1) = Cx(1, 1);
    vaux(1, 2) = Cy(1, 1);
    vaux(1, 3) = Cz(1, 1);
    Matrix ChebApp = (f1 * tau) - f2 + vaux;

    return ChebApp;
}