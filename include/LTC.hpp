// $Header$
//--------------------------------------------------
// LTC
//--------------------------------------------------
// Proyecto_TT1: Local Tangent Coordinates Transformation
//
// Copyright (c) 2025, Meysam Mahooti
//
// Created: 2025/05/07
//
/** @file LTC.hpp
 *  @brief Transformation from Greenwich meridian system to local tangent coordinates
 *
 *  @author Noel Del Rio Gonzalez
 *  @bug No known bugs
 */

#ifndef _LTC_
#define _LTC_

#include "matrix.hpp"

Matrix LTC(double lon, double lat);

#endif