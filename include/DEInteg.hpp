// $Source$
//--------------------------------------------------
// DEInteg
//--------------------------------------------------
// Proyecto_TT1: Numerical Integration for ODEs
//
// Copyright (c) 2025, Meysam Mahooti
//
// Created: 2025/05/23
//
/** @file DEInteg.hpp
 *  @brief Header for numerical integration of ordinary differential equations
 *
 *  @author NOEL DEL RIO GONZALEZ
 *  @bug No known bugs
 */

#ifndef DEINTEG_HPP
#define DEINTEG_HPP

#include "..\include\AccelPointMass.hpp"
#include "..\include\Cheb3D.hpp"
#include "..\include\EccAnom.hpp"
#include "..\include\Frac.hpp"
#include "..\include\MeanObliquity.hpp"
#include "..\include\Mjday.hpp"
#include "..\include\Mjday_TDB.hpp"
#include "..\include\Position.hpp"
#include "..\include\R_x.hpp"
#include "..\include\R_y.hpp"
#include "..\include\R_z.hpp"
#include "..\include\sign_.hpp"
#include "..\include\timediff.hpp"
#include "..\include\AzElPa.hpp"
#include "..\include\IERS.hpp"
#include "..\include\Legendre.hpp"
#include "..\include\NutAngles.hpp"
#include "..\include\TimeUpdate.hpp"
#include "..\include\LTC.hpp"
#include "..\include\NutMatrix.hpp"
#include "..\include\PoleMatrix.hpp"
#include "..\include\PrecMatrix.hpp"
#include "..\include\gmst.hpp"
#include "..\include\AccelHarmonic.hpp"
#include "..\include\EqnEquinox.hpp"
#include "..\include\JPL_Eph_DE430.hpp"
#include "..\include\NutMatrix.hpp"
#include "..\include\PrecMatrix.hpp"
#include "..\include\GMST.hpp"
#include "..\include\gast.hpp"
#include "..\include\MeasUpdate.hpp"
#include "..\include\G_AccelHarmonic.hpp"
#include "..\include\GHAMatrix.hpp"
#include "..\include\Accel.hpp"
#include "..\include\VarEqn.hpp"



#include "Matrix.hpp"

/**
 * @brief Numerical integration for ordinary differential equations
 * @param f Function defining the ODE (dy/dt = f(t,y))
 * @param t Initial time (seconds)
 * @param tout Desired output time (seconds)
 * @param relerr Relative error tolerance
 * @param abserr Absolute error tolerance
 * @param n_eqn Number of equations
 * @param y Initial state vector (n_eqn x 1)
 * @return y Updated state vector at tout (n_eqn x 1)
 * @throw std::invalid_argument If parameters are invalid
 */
Matrix DEInteg(Matrix (*f)(double, Matrix&), double t, double tout, double relerr, double abserr, int n_eqn, Matrix& y);


#endif // DEINTEG_HPP