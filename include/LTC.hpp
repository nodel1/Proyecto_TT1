// $Header$
//--------------------------------------------------
// LTC
//--------------------------------------------------
// Proyecto_TT1: Local Tangent Coordinates Transformation Header
//
// Copyright (c) 2025, Meysam Mahooti
//
// Created: 2025/05/08
//
/** @file LTC.hpp
 *  @brief Header for transformation from Greenwich meridian system to local tangent coordinates
 *
 *  @author Noel Del Rio Gonzalez
 *  @bug No known bugs
 */

#ifndef LTC_HPP
#define LTC_HPP

#include "matrix.hpp"

Matrix LTC(double lon, double lat);

#endif