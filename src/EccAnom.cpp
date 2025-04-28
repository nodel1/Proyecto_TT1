// $Source$
//--------------------------------------------------
// EccAnom
//--------------------------------------------------
// Proyecto_TT1: Eccentric Anomaly Calculation Implementation
//
// Copyright (c) 2025, Meysam Mahooti
//
// Created: 2025/04/23
//
/** @file ecc_anom.cpp
 *  @brief Eccentric anomaly computation for elliptic orbits
 *
 *  @author Noel Del Rio Gonzalez
 *  @bug No known bugs
 */

#include "..\include\EccAnom.hpp"

using namespace std;

double EccAnom(double M, double e) {
    const int maxit = 15;
    int i = 1;
    double E, f;

    // Starting value
    M = fmod(M, 2.0 * M_PI);
    E = (e < 0.8) ? M : M_PI;

    // Initial calculation
    f = E - e * sin(E) - M;
    E = E - f / (1.0 - e * cos(E));

    // Iteration
    while (abs(f) > 1e2 * numeric_limits<double>::epsilon()) {
        f = E - e * sin(E) - M;
        E = E - f / (1.0 - e * cos(E));
        if (++i > maxit) {
            cout << "EccAnom: convergence problems after " << maxit << " iterations\n";
            exit(EXIT_FAILURE);
        }
    }

    return E;
}