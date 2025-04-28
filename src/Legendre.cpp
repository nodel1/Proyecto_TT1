// $Source$
//--------------------------------------------------
// Legendre
//--------------------------------------------------
// Proyecto_TT1: Associated Legendre Functions Implementation
//
// Copyright (c) 2025, Meysam Mahooti
//
// Created: 2025/04/23
//
/** @file Legendre.cpp
 *  @brief Implementation of normalized associated Legendre functions and their derivatives
 *
 *  @author Noel Del Rio Gonzalez
 *  @bug No known bugs
 */

#include "..\include\Legendre.hpp"

using namespace std;

pair<Matrix, Matrix> Legendre(int n, int m, double fi) {
    Matrix pnm = zeros(n + 1, m + 1);   // Initialize (n+1) x (m+1) matrix with zeros
    Matrix dpnm = zeros(n + 1, m + 1);  // Initialize (n+1) x (m+1) matrix with zeros

    // Initial values
    pnm(0, 0) = 1.0;
    dpnm(0, 0) = 0.0;
    if (n >= 1 && m >= 1) {
        pnm(1, 1) = sqrt(3.0) * cos(fi);
        dpnm(1, 1) = -sqrt(3.0) * sin(fi);
    }

    // Diagonal coefficients
    for (int i = 2; i <= n; i++) {
        pnm(i, i) = sqrt((2.0 * i + 1.0) / (2.0 * i)) * cos(fi) * pnm(i - 1, i - 1);
    }
    for (int i = 2; i <= n; i++) {
        dpnm(i, i) = sqrt((2.0 * i + 1.0) / (2.0 * i)) * (
            cos(fi) * dpnm(i - 1, i - 1) - sin(fi) * pnm(i - 1, i - 1)
        );
    }

    // Horizontal first step coefficients
    for (int i = 1; i <= n; i++) {
        pnm(i, i - 1) = sqrt(2.0 * i + 1.0) * sin(fi) * pnm(i - 1, i - 1);
    }
    for (int i = 1; i <= n; i++) {
        dpnm(i, i - 1) = sqrt(2.0 * i + 1.0) * (
            cos(fi) * pnm(i - 1, i - 1) + sin(fi) * dpnm(i - 1, i - 1)
        );
    }

    // Horizontal second step coefficients
    int j = 0;
    int k = 2;
    while (true) {
        for (int i = k; i <= n; i++) {
            pnm(i, j) = sqrt((2.0 * i + 1.0) / ((i - j) * (i + j))) * (
                sqrt(2.0 * i - 1.0) * sin(fi) * pnm(i - 1, j) -
                sqrt(((i + j - 1.0) * (i - j - 1.0)) / (2.0 * i - 3.0)) * pnm(i - 2, j)
            );
        }
        j++;
        k++;
        if (j > m) {
            break;
        }
    }

    j = 0;
    k = 2;
    while (true) {
        for (int i = k; i <= n; i++) {
            dpnm(i, j) = sqrt((2.0 * i + 1.0) / ((i - j) * (i + j))) * (
                sqrt(2.0 * i - 1.0) * sin(fi) * dpnm(i - 1, j) +
                sqrt(2.0 * i - 1.0) * cos(fi) * pnm(i - 1, j) -
                sqrt(((i + j - 1.0) * (i - j - 1.0)) / (2.0 * i - 3.0)) * dpnm(i - 2, j)
            );
        }
        j++;
        k++;
        if (j > m) {
            break;
        }
    }

    return {pnm, dpnm};
}