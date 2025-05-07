// $Header$
//--------------------------------------------------
// EqnEquinox
//--------------------------------------------------
// Proyecto_TT1: Equation of the Equinoxes Header
//
// Copyright (c) 2025, Meysam Mahooti
//
// Created: 2025/05/07
//
/** @file EqnEquinox.hpp
 *  @brief Computes the equation of the equinoxes
 *
 *  @author Noel Del Rio Gonzalez
 *  @bug No known bugs
 */

#ifndef EQNEQUINOX_HPP
#define EQNEQUINOX_HPP

#include "NutAngles.hpp"
#include "MeanObliquity.hpp"

//--------------------------------------------------
// EqnEquinox (double Mjd_TT)
//--------------------------------------------------
/**
 * @brief Computes the equation of the equinoxes
 *
 * @param[in] Mjd_TT Modified Julian Date (Terrestrial Time)
 * @return double Equation of the equinoxes (EqE) [rad]
 */
double EqnEquinox(double Mjd_TT);

#endif // EQNEQUINOX_HPP