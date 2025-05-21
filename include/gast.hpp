// $Source$
//--------------------------------------------------
// GAST
//--------------------------------------------------
// Proyecto_TT1: Greenwich Apparent Sidereal Time Declaration
//
// Copyright (c) 2025, Meysam Mahooti
//
// Created: 2025/05/08
//
/** @file GAST.hpp
 *  @brief Declaration of Greenwich Apparent Sidereal Time calculation
 *
 *  @author Noel Del Rio Gonzalez
 *  @bug No known bugs
 */

#ifndef GAST_HPP
#define GAST_HPP

/**
 * @brief Calculate Greenwich Apparent Sidereal Time
 * @param Mjd_UT1 Modified Julian Date (UT1)
 * @return GAST in radians [0, 2pi)
 */
double GAST(double Mjd_UT1);

#endif // GAST_HPP