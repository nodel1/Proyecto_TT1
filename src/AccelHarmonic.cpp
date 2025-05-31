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
    std::cout << "Starting AccelHarmonic computation..." << std::endl;
    
    // Constantes
    const double r_ref = 6378.1363e3;   // Earth's radius [m]; GGM03S
    const double gm = 398600.4415e9;    // [m^3/s^2]; GGM03S
    std::cout << "Constants loaded: r_ref = " << r_ref << " m, gm = " << gm << " m^3/s^2" << std::endl;

    // Body-fixed position
    std::cout << "Calculating body-fixed position..." << std::endl;
    Matrix r_bf = E * r;
    std::cout << "Body-fixed position calculated:\n" << r_bf << std::endl;

    // Auxiliary quantities
    std::cout << "Calculating auxiliary quantities..." << std::endl;
    double d = norm(r_bf);
    double latgc = std::asin(r_bf(3, 1) / d);
    double lon = std::atan2(r_bf(2, 1), r_bf(1, 1));
    std::cout << "Auxiliary quantities: d = " << d << " m, latgc = " << latgc 
              << " rad, lon = " << lon << " rad" << std::endl;

    // Calcular funciones de Legendre
    std::cout << "Calculating Legendre functions..." << std::endl;
    Matrix pnm(n_max + 1, m_max + 1);
    Matrix dpnm(n_max + 1, m_max + 1);
	
	
	    std::cout << "Printing pnm matrix (" << n_max + 1 << " x " << m_max + 1 << "):" << std::endl;
    for (int i = 0; i <= n_max; i++) {
        for (int j = 0; j <= m_max; j++) {
            std::cout << "pnm(" << i << "," << j << ") = " << pnm(i+1, j+1) << "  ";
        }
        std::cout << std::endl;
    }
    std::cout << "Printing dpnm matrix (" << n_max + 1 << " x " << m_max + 1 << "):" << std::endl;
    for (int i = 0; i <= n_max; i++) {
        for (int j = 0; j <= m_max; j++) {
            std::cout << "dpnm(" << i << "," << j << ") = " << dpnm(i+1, j+1) << "  ";
        }
        std::cout << std::endl;
    }
	
	
	
	
    Legendre(n_max, m_max, latgc, pnm, dpnm);
    std::cout << "Legendre functions calculated (pnm and dpnm matrices)" << std::endl;
	
	
	
	    std::cout << "Printing pnm matrix after Legendre (" << n_max + 1 << " x " << m_max + 1 << "):" << std::endl;
    for (int i = 0; i <= n_max; i++) {
        for (int j = 0; j <= m_max; j++) {
            std::cout << "pnm(" << i << "," << j << ") = " << pnm(i+1, j+1) << "  ";
        }
        std::cout << std::endl;
    }
    std::cout << "Printing dpnm matrix after Legendre (" << n_max + 1 << " x " << m_max + 1 << "):" << std::endl;
    for (int i = 0; i <= n_max; i++) {
        for (int j = 0; j <= m_max; j++) {
            std::cout << "dpnm(" << i << "," << j << ") = " << dpnm(i+1, j+1) << "  ";
        }
        std::cout << std::endl;
    }
	
	
	

    // Inicializar acumuladores
    std::cout << "Initializing accumulators for potential derivatives..." << std::endl;
    double dUdr = 0.0;
    double dUdlatgc = 0.0;
    double dUdlon = 0.0;
    double q1 = 0.0, q2 = 0.0, q3 = 0.0;
	
	
	if (n_max > 10){
	std::cout << "pnm(4,4) = " << pnm(4,4) << ", dpnm(4,4) = " << dpnm(4,4) << std::endl;
	std::cout << "pnm(4,4) = " << pnm(4,4) << ", dpnm(4,4) = " << dpnm(4,4) << std::endl;
	std::cout << "pnm(4,4) = " << pnm(4,4) << ", dpnm(4,4) = " << dpnm(4,4) << std::endl;
	std::cout << "pnm(4,4) = " << pnm(4,4) << ", dpnm(4,4) = " << dpnm(4,4) << std::endl;
	std::cout << "pnm(4,4) = " << pnm(4,4) << ", dpnm(4,4) = " << dpnm(4,4) << std::endl;
	std::cout << "pnm(4,4) = " << pnm(4,4) << ", dpnm(4,4) = " << dpnm(4,4) << std::endl;
	std::cout << "pnm(4,4) = " << pnm(4,4) << ", dpnm(4,4) = " << dpnm(4,4) << std::endl;
	std::cout << "pnm(4,4) = " << pnm(4,4) << ", dpnm(4,4) = " << dpnm(4,4) << std::endl;
	std::cout << "pnm(4,4) = " << pnm(4,4) << ", dpnm(4,4) = " << dpnm(4,4) << std::endl;
	std::cout << "pnm(4,4) = " << pnm(4,4) << ", dpnm(4,4) = " << dpnm(4,4) << std::endl;
	std::cout << "pnm(4,4) = " << pnm(4,4) << ", dpnm(4,4) = " << dpnm(4,4) << std::endl;
	std::cout << "pnm(4,4) = " << pnm(4,4) << ", dpnm(4,4) = " << dpnm(4,4) << std::endl;
	std::cout << "pnm(4,4) = " << pnm(4,4) << ", dpnm(4,4) = " << dpnm(4,4) << std::endl;
	std::cout << "pnm(4,4) = " << pnm(4,4) << ", dpnm(4,4) = " << dpnm(4,4) << std::endl;
	std::cout << "pnm(4,4) = " << pnm(4,4) << ", dpnm(4,4) = " << dpnm(4,4) << std::endl;
	std::cout << "pnm(4,4) = " << pnm(4,4) << ", dpnm(4,4) = " << dpnm(4,4) << std::endl;
	std::cout << "pnm(0,3) = " << pnm(0,3) << ", dpnm(0,3) = " << dpnm(0,3) << std::endl;
	std::cout << "pnm(0,4) = " << pnm(0,4) << ", dpnm(0,4) = " << dpnm(0,4) << std::endl;
	}
	
	

    // Bucle para calcular dUdr, dUdlatgc, dUdlon
    std::cout << "Starting harmonic summation (n_max = " << n_max << ", m_max = " << m_max << ")..." << std::endl;
    for (int n = 0; n <= n_max; n++) {
        double b1 = (-gm / (d * d)) * std::pow(r_ref / d, n) * (n + 1);
        double b2 = (gm / d) * std::pow(r_ref / d, n);
        double b3 = (gm / d) * std::pow(r_ref / d, n);
        q1 = 0.0; q2 = 0.0; q3 = 0.0;
        std::cout << "AQUI LLEGAS "<< std::endl;
        for (int m = 0; m <= m_max; m++) {
				std::cout << "pnm(" << n << "," << m << ") = " << pnm(n+1,m+1) << ", dpnm(" << n << "," << m << ") = " << dpnm(n+1,m+1) << std::endl;
			        std::cout << "AQUI LLEGASss "<< std::endl;	
            q1 += pnm(n + 1, m + 1) * (Cnm(n + 1, m + 1) * std::cos(m * lon) +
                                     Snm(n + 1, m + 1) * std::sin(m * lon));
			
			std::cout << "AQUI LLEGAS2 ademas con una m,n con valor" << std::endl;
			
            q2 += dpnm(n + 1, m + 1) * (Cnm(n + 1, m + 1) * std::cos(m * lon) +
                                      Snm(n + 1, m + 1) * std::sin(m * lon));
            q3 += m * pnm(n + 1, m + 1) * (Snm(n + 1, m + 1) * std::cos(m * lon) -
                                         Cnm(n + 1, m + 1) * std::sin(m * lon));
        }

        dUdr += q1 * b1;
        dUdlatgc += q2 * b2;
        dUdlon += q3 * b3;
        
        if (n % 5 == 0) {  // Progress report every 5 degrees
            std::cout << "  Completed degree " << n << "/" << n_max 
                      << " (dUdr = " << dUdr << ", dUdlatgc = " << dUdlatgc 
                      << ", dUdlon = " << dUdlon << ")" << std::endl;
        }
    }
    std::cout << "Harmonic summation completed" << std::endl;
    std::cout << "Final potential derivatives: dUdr = " << dUdr 
              << ", dUdlatgc = " << dUdlatgc << ", dUdlon = " << dUdlon << std::endl;

    // Body-fixed acceleration
    std::cout << "Calculating body-fixed acceleration..." << std::endl;
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
    std::cout << "Body-fixed acceleration calculated:\n" << a_bf << std::endl;

    // Inertial acceleration
    std::cout << "Transforming to inertial acceleration..." << std::endl;
    Matrix E_t = transpose(E);
    Matrix a = E_t * a_bf;
    std::cout << "Inertial acceleration calculated:\n" << a << std::endl;

    std::cout << "AccelHarmonic computation completed successfully" << std::endl;
    return a;
}