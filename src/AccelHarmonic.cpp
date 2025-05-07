// $Source$
//--------------------------------------------------
// AccelHarmonic
//--------------------------------------------------
// Proyecto_TT1: Harmonic Gravity Field Acceleration Implementation
//
// Copyright (c) 2025, Meysam Mahooti
//
// Created: 2025/05/07
//
/** @file AccelHarmonic.cpp
 *  @brief Computes the acceleration due to the harmonic gravity field of the central body
 *
 *  @author Noel Del Rio Gonzalez
 *  @bug No known bugs
 */

#include "..\include\AccelHarmonic.hpp"
#include "..\include\SAT_Const.hpp"
#include "..\include\global.hpp"
#include <cmath>


Matrix AccelHarmonic(Matrix& r, Matrix& E, int n_max, int m_max) {
    // Constantes
    const double r_ref = 6378.1363e3;   // Earth's radius [m]; GGM03S
    const double gm = 398600.4415e9;    // [m^3/s^2]; GGM03S

    // Body-fixed position
    Matrix r_bf = E * r;

    // Auxiliary quantities
    double d = norm(r_bf);
    double latgc = std::asin(r_bf(3, 1) / d);
    double lon = std::atan2(r_bf(2, 1), r_bf(1, 1));

    // Calcular funciones de Legendre
    Matrix pnm(n_max + 1, m_max + 1);
    Matrix dpnm(n_max + 1, m_max + 1);        //REVISAR ESTO 
    Legendre(n_max, m_max, latgc, pnm, dpnm);

    // Inicializar acumuladores
    double dUdr = 0.0;
    double dUdlatgc = 0.0;
    double dUdlon = 0.0;
    double q1 = 0.0, q2 = 0.0, q3 = 0.0;

    // Bucle para calcular dUdr, dUdlatgc, dUdlon
    for (int n = 0; n <= n_max; n++) {
        double b1 = (-gm / (d * d)) * std::pow(r_ref / d, n) * (n + 1);
        double b2 = (gm / d) * std::pow(r_ref / d, n);
        double b3 = (gm / d) * std::pow(r_ref / d, n);
        q1 = 0.0; q2 = 0.0; q3 = 0.0;

        for (int m = 0; m <= m_max; m++) {
            q1 += pnm(n + 1, m + 1) * (Cnm(n + 1, m + 1) * std::cos(m * lon) +
                                       Snm(n + 1, m + 1) * std::sin(m * lon));
            q2 += dpnm(n + 1, m + 1) * (Cnm(n + 1, m + 1) * std::cos(m * lon) +
                                        Snm(n + 1, m + 1) * std::sin(m * lon));
            q3 += m * pnm(n + 1, m + 1) * (Snm(n + 1, m + 1) * std::cos(m * lon) -
                                           Cnm(n + 1, m + 1) * std::sin(m * lon));
        }

        dUdr += q1 * b1;
        dUdlatgc += q2 * b2;
        dUdlon += q3 * b3;
    }

    // Body-fixed acceleration
    double r2xy = r_bf(1, 1) * r_bf(1, 1) + r_bf(2, 1) * r_bf(2, 1);
    double ax = (1.0 / d * dUdr - r_bf(3, 1) / (d * d * std::sqrt(r2xy)) * dUdlatgc) * r_bf(1, 1) -
                (1.0 / r2xy * dUdlon) * r_bf(2, 1);
    double ay = (1.0 / d * dUdr - r_bf(3, 1) / (d * d * std::sqrt(r2xy)) * dUdlatgc) * r_bf(2, 1) +
                (1.0 / r2xy * dUdlon) * r_bf(1, 1);
    double az = 1.0 / d * dUdr * r_bf(3, 1) + std::sqrt(r2xy) / (d * d) * dUdlatgc;

    Matrix a_bf(3, 1);
    a_bf(1, 1) = ax;
    a_bf(2, 1) = ay;
    a_bf(3, 1) = az;

    // Inertial acceleration
    Matrix E_t = transpose(E);
    Matrix a = E_t * a_bf;

    return a;
}