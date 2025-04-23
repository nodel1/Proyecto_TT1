// $Header$
//--------------------------------------------------
// MeanObliquity
//--------------------------------------------------
// Proyecto_TT1: Mean Obliquity Calculation Header
//
// Copyright (c) 2025, Meysam Mahooti
//
// Created: 2025/04/23
//
/** @file mean_obliquity.hpp
 *  @brief Mean obliquity of the ecliptic computation
 *
 *  @author Noel Del Rio Gonzalez
 *  @bug No known bugs
 */

#ifndef MEAN_OBLIQUITY_HPP
#define MEAN_OBLIQUITY_HPP

//--------------------------------------------------
// MeanObliquity (double Mjd_TT)
//--------------------------------------------------
/**
 * @brief Computes the mean obliquity of the ecliptic
 * 
 * @param[in] Mjd_TT Modified Julian Date (Terrestrial Time)
 * @return double Mean obliquity of the ecliptic [rad]
 * 
 * @note Uses IAU 2006 precession model
 * @note Requires global const struct with constants:
 *       - MJD_J2000
 *       - Rad
 */
double MeanObliquity(double Mjd_TT);

#endif // MEAN_OBLIQUITY_HPP