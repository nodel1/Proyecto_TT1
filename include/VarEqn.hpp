// $Source$
//--------------------------------------------------
// VarEqn
//--------------------------------------------------
// Proyecto_TT1: Variational Equations Computation
//
// Copyright (c) 2025, Meysam Mahooti
//
// Created: 2025/05/23
//
/** @file VarEqn.hpp
 *  @brief Computes variational equations for orbit propagation
 *
 *  @author Noel Del Rio
 *  @bug No known bugs
 */

#ifndef VAREQN_HPP
#define VAREQN_HPP

#include "Matrix.hpp"

/**
 * @brief Computes the derivative of the state vector and transition matrix
 * @param x Time since epoch in seconds
 * @param yPhi Combined state vector and state transition matrix (42x1)
 * @return Derivative vector (42x1)
 */
Matrix VarEqn(double x, Matrix yPhi);

#endif // VAREQN_HPP
