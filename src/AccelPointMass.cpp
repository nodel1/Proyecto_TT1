// $Source$
//--------------------------------------------------
// AccelPointMass
//--------------------------------------------------
// Proyecto_TT1: Point Mass Acceleration Implementation
//
// Copyright (c) 2025, Meysam Mahooti
//
// Created: 2025/04/23
//
/** @file accel_point_mass.cpp
 *  @brief Implementation of point mass perturbation acceleration computation
 *
 *  @author Noel Del Rio Gonzalez
 *  @bug No known bugs
 */

#include "..\include\AccelPointMass.hpp"
#include <cmath>

Matrix& AccelPointMass(Matrix& r, Matrix& s, double GM) {
    // Relative position vector of satellite w.r.t. point mass
	    std::cout << "antes de la resta de matrices" << std::endl;
    Matrix& d = r - s;
    	    std::cout << "despues" << std::endl;
    // Acceleration
    return (d/pow(norm(d),3) + s/pow(norm(s),3)) * (-GM);
	    	    std::cout << "salida de accelpointmass" << std::endl;
}