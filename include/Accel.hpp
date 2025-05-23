// $Source$
//--------------------------------------------------
// Accel
//--------------------------------------------------
// Proyecto_TT1: Satellite Acceleration Computation
//
// Copyright (c) 2025, Meysam Mahooti
//
// Created: 2025/05/23
//
/** @file Accel.hpp
 *  @brief Computes the acceleration of an Earth-orbiting satellite
 *
 *  @author Noel Del Rio
 *  @bug No known bugs
 */

#ifndef ACCEL_HPP
#define ACCEL_HPP

#include "Matrix.hpp"

/**
 * @brief Computes the acceleration of an Earth-orbiting satellite
 * @param x Terrestrial Time (seconds since Modified Julian Date)
 * @param Y Satellite state vector in the ICRF/EME2000 system (6x1)
 * @return dY Acceleration vector (6x1) in the ICRF/EME2000 system
 */
Matrix Accel(double x, Matrix Y);

#endif // ACCEL_HPP