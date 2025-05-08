// $Header$
//--------------------------------------------------
// JPL_Eph_DE430
//--------------------------------------------------
// Proyecto_TT1: JPL DE430 Ephemerides Header
//
// Copyright (c) 2025, Meysam Mahooti
//
// Created: 2025/05/07
//
/** @file JPL_Eph_DE430.hpp
 *  @brief Computes the sun, moon, and nine major planets' equatorial position using JPL DE430 Ephemerides
 *
 *  @author Noel Del Rio Gonzalez
 *  @bug No known bugs
 */

#ifndef JPL_EPH_DE430_HPP
#define JPL_EPH_DE430_HPP

#include "global.hpp"
#include "matrix.hpp"
#include "Cheb3D.hpp"

//--------------------------------------------------
// JPL_Eph_DE430 (double Mjd_TDB, Matrix& r_Mercury, ...)
//--------------------------------------------------
/**
 * @brief Computes the equatorial positions of the sun, moon, and nine major planets using JPL DE430 Ephemerides
 *
 * @param[in] Mjd_TDB Modified Julian Date (TDB)
 * @param[out] r_Mercury Geocentric equatorial position of Mercury [m] (3x1)
 * @param[out] r_Venus Geocentric equatorial position of Venus [m] (3x1)
 * @param[out] r_Earth Solar system barycenter position of Earth [m] (3x1)
 * @param[out] r_Mars Geocentric equatorial position of Mars [m] (3x1)
 * @param[out] r_Jupiter Geocentric equatorial position of Jupiter [m] (3x1)
 * @param[out] r_Saturn Geocentric equatorial position of Saturn [m] (3x1)
 * @param[out] r_Uranus Geocentric equatorial position of Uranus [m] (3x1)
 * @param[out] r_Neptune Geocentric equatorial position of Neptune [m] (3x1)
 * @param[out] r_Pluto Geocentric equatorial position of Pluto [m] (3x1)
 * @param[out] r_Moon Geocentric equatorial position of Moon [m] (3x1)
 * @param[out] r_Sun Geocentric equatorial position of Sun [m] (3x1)
 *
 */
void JPL_Eph_DE430(double Mjd_TDB, Matrix& r_Mercury, Matrix& r_Venus, Matrix& r_Earth,
                   Matrix& r_Mars, Matrix& r_Jupiter, Matrix& r_Saturn, Matrix& r_Uranus,
                   Matrix& r_Neptune, Matrix& r_Pluto, Matrix& r_Moon, Matrix& r_Sun);

#endif // JPL_EPH_DE430_HPP
