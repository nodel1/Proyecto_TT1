// $Source$
//--------------------------------------------------
// TimeUpdate
//--------------------------------------------------
// Proyecto_TT1: Covariance Matrix Update Implementation
//
// Copyright (c) 2025, Meysam Mahooti
//
// Created: 2025/04/23
//
/** @file TimeUpdate.cpp
 *  @brief Implementation of covariance matrix update
 *
 *  @author Noel Del Rio Gonzalez
 *  @bug No known bugs
 */

#include "TimeUpdate.hpp"

using namespace std;

Matrix TimeUpdate(const Matrix& P, const Matrix& Phi, const Matrix& Qdt) {
    Matrix P_new = Phi * P * transpose(Phi) + Qdt;
    return P_new;
}

Matrix TimeUpdate(const Matrix& P, const Matrix& Phi) {
    Matrix Qdt = zeros(P.n_row, P.n_row); // Zero process noise matrix
    Matrix P_new = Phi * P * transpose(Phi) + Qdt;
    return P_new;
}