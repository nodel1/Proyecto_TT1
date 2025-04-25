// $Source$
//--------------------------------------------------
// sign_
//--------------------------------------------------
// Proyecto_TT1: Sign Function Implementation
//
// Copyright (c) 2025, Meysam Mahooti
//
// Created: 2025/04/23
//
/** @file sign_.cpp
 *  @brief Implementation of the sign function
 *
 *  @author Noel Del Rio Gonzalez
 *  @bug No known bugs
 */

#include "sign_.hpp"

using namespace std;

double sign_(double a, double b) {
    if (b >= 0.0) {
        return abs(a);
    } else {
        return -abs(a);
    }
}