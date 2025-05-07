// $Header$
//--------------------------------------------------
// AccelHarmonic
//--------------------------------------------------
// Proyecto_TT1: Harmonic Gravity Field Acceleration Header
//
// Copyright (c) 2025, Meysam Mahooti
//
// Created: 2025/05/07
//
/** @file AccelHarmonic.hpp
 *  @brief Computes the acceleration due to the harmonic gravity field of the central body
 *
 *  @author Noel Del Rio Gonzalez
 *  @bug No known bugs
 */

#ifndef ACCELHARMONIC_HPP
#define ACCELHARMONIC_HPP

#include "matrix.hpp"
#include "Legendre.hpp"
#include "global.hpp"
#include <cmath>

//--------------------------------------------------
// AccelHarmonic (Matrix& r, Matrix& E, int n_max, int m_max)
//--------------------------------------------------
/**
 * @brief Computes the acceleration due to the harmonic gravity field of the central body
 *
 * @param[in] r Satellite position vector in the inertial system (3x1)
 * @param[in] E Transformation matrix to body-fixed system (3x3)
 * @param[in] n_max Maximum degree
 * @param[in] m_max Maximum order (m_max<=n_max; m_max=0 for zonals, only)
 * @return Matrix Acceleration (a=d^2r/dt^2) (3x1)
 *
 */
Matrix AccelHarmonic(Matrix& r, Matrix& E, int n_max, int m_max);

#endif // ACCELHARMONIC_HPP