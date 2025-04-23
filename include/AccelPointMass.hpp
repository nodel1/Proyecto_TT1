// $Header$
//--------------------------------------------------
// AccelPointMass
//--------------------------------------------------
// Proyecto_TT1: Point Mass Acceleration Header
//
// Copyright (c) 2025, Meysam Mahooti
//
// Created: 2025/04/23
//
/** @file accel_point_mass.hpp
 *  @brief Header for point mass perturbation acceleration computation
 *
 *  @author Noel Del Rio Gonzalez
 *  @bug No known bugs
 */

#ifndef ACCEL_POINT_MASS_HPP
#define ACCEL_POINT_MASS_HPP

#include "matrix.hpp"

//--------------------------------------------------
// AccelPointMass (Matrix &r, Matrix &s, double GM)
//--------------------------------------------------
/**
 * @brief Computes the perturbational acceleration due to a point mass
 * 
 * @param[in] r Satellite position vector (3x1 matrix)
 * @param[in] s Point mass position vector (3x1 matrix)
 * @param[in] GM Gravitational coefficient of point mass
 * @return Matrix& Acceleration vector (3x1 matrix)
 * 
 * @note Both input vectors must be 3x1 matrices
 */
Matrix& AccelPointMass(Matrix& r, Matrix& s, double GM);

#endif // ACCEL_POINT_MASS_HPP