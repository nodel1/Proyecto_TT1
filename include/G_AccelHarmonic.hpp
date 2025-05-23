// $Source$
//--------------------------------------------------
// G_AccelHarmonic
//--------------------------------------------------
// Proyecto_TT1: Gradient of Harmonic Gravity Field Implementation
//
// Copyright (c) 2025, Meysam Mahooti
//
// Created: 2025/05/23
//
/** @file G_AccelHarmonic.hpp
 *  @brief Computes the gradient of the Earth's harmonic gravity field
 *
 *  @author Noel Del Rio Gonzalez
 *  @bug No known bugs
 */

#ifndef G_ACCELHARMONIC_HPP
#define G_ACCELHARMONIC_HPP

#include "Matrix.hpp"

/**
 * @brief Computes the gradient of the Earth's harmonic gravity field
 *
 * @param r Satellite position vector in the true-of-date system (3x1 matrix)
 * @param U Transformation matrix to body-fixed system (3x3 matrix)
 * @param n_max Gravity model degree
 * @param m_max Gravity model order (m_max <= n_max)
 * @return G Gradient (G = da/dr) in the true-of-date system (3x3 matrix)
 */
Matrix G_AccelHarmonic(Matrix& r, Matrix& U, int n_max, int m_max);

#endif // G_ACCELHARMONIC_HPP