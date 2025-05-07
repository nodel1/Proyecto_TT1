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

void Legendre(int n, int m, double fi, Matrix& pnm, Matrix& dpnm) {
    pnm = zeros(n+1, m+1);
    dpnm = zeros(n+1, m+1);

    // Initial values
    pnm(1,1) = 1.0;
    dpnm(1,1) = 0.0;
    if (n >= 1 && m >= 1) {
        pnm(2,2) = std::sqrt(3.0) * std::cos(fi);
        dpnm(2,2) = -std::sqrt(3.0) * std::sin(fi);
    }

    // Diagonal coefficients
    for (int i = 2; i <= n; i++) {
        pnm(i+1,i+1) = std::sqrt((2.0*i+1)/(2.0*i)) * std::cos(fi) * pnm(i,i);
    }
    for (int i = 2; i <= n; i++) {
        dpnm(i+1,i+1) = std::sqrt((2.0*i+1)/(2.0*i)) * (
            std::cos(fi) * dpnm(i,i) - std::sin(fi) * pnm(i,i)
        );
    }

    // Horizontal first step coefficients
    for (int i = 1; i <= n; i++) {
        pnm(i+1,i) = std::sqrt(2.0*i+1) * std::sin(fi) * pnm(i,i);
    }
    for (int i = 1; i <= n; i++) {
        dpnm(i+1,i) = std::sqrt(2.0*i+1) * (
            std::cos(fi) * pnm(i,i) + std::sin(fi) * dpnm(i,i)
        );
    }

    // Horizontal second step coefficients
    int j = 0;
    int k = 2;
    while (true) {
        for (int i = k; i <= n; i++) {
            pnm(i+1,j+1) = std::sqrt((2.0*i+1)/((i-j)*(i+j))) * (
                std::sqrt(2.0*i-1) * std::sin(fi) * pnm(i,j+1) -
                std::sqrt(((i+j-1.0)*(i-j-1.0))/(2.0*i-3)) * pnm(i-1,j+1)
            );
        }
        j++;
        k++;
        if (j > m) break;
    }
    j = 0;
    k = 2;
    while (true) {
        for (int i = k; i <= n; i++) {
            dpnm(i+1,j+1) = std::sqrt((2.0*i+1)/((i-j)*(i+j))) * (
                std::sqrt(2.0*i-1) * std::sin(fi) * dpnm(i,j+1) +
                std::sqrt(2.0*i-1) * std::cos(fi) * pnm(i,j+1) -
                std::sqrt(((i+j-1.0)*(i-j-1.0))/(2.0*i-3)) * dpnm(i-1,j+1)
            );
        }
        j++;
        k++;
        if (j > m) break;
    }
}