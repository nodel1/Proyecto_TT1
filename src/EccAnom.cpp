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
/** @file EccAnom.cpp
 *  @brief Eccentric anomaly computation for elliptic orbits
 *
 *  @author Noel Del Rio Gonzalez
 *  @bug No known bugs
 */

#include <cmath>
#include <numbers>
#include <stdexcept>
#include <limits>
#include "..\include\EccAnom.hpp"

double EccAnom(double M, double e) {
    const int maxit = 15;
    int i = 1;
    double E, f;

    // Starting value
    M = std::fmod(M, 2.0 * std::numbers::pi);
    E = (e < 0.8) ? M : std::numbers::pi;

    // Initial calculation
    f = E - e * std::sin(E) - M;
    E = E - f / (1.0 - e * std::cos(E));

    // Iteration
    while (std::abs(f) > 1e2 * std::numeric_limits<double>::epsilon()) {
        f = E - e * std::sin(E) - M;
        E = E - f / (1.0 - e * std::cos(E));
        if (++i > maxit) {
            throw std::runtime_error("EccAnom: convergence problems after " + std::to_string(maxit) + " iterations");
        }
    }

    return E;
}