// $Header$
//--------------------------------------------------
// Cheb3D
//--------------------------------------------------
// Proyecto_TT1: Chebyshev Approximation Header
//
// Copyright (c) 2025, Meysam Mahooti
//
// Created: 2025/04/23
//
/** @file cheb3d.hpp
 *  @brief 3D Chebyshev approximation implementation
 *
 *  @author Noel Del Rio Gonzalez
 *  @bug No known bugs
 */

#ifndef CHEB3D_HPP
#define CHEB3D_HPP

#include "matrix.hpp"

//--------------------------------------------------
// Cheb3D (double t, int N, double Ta, double Tb, 
//         Matrix &Cx, Matrix &Cy, Matrix &Cz)
//--------------------------------------------------
/**
 * @brief Computes 3D Chebyshev approximation
 * 
 * @param[in] t Evaluation time
 * @param[in] N Number of coefficients
 * @param[in] Ta Begin interval
 * @param[in] Tb End interval
 * @param[in] Cx Coefficients for x-coordinate (Nx1 matrix)
 * @param[in] Cy Coefficients for y-coordinate (Nx1 matrix)
 * @param[in] Cz Coefficients for z-coordinate (Nx1 matrix)
 * @return Matrix& Approximated vector (1x3 matrix)
 * 
 */
 
     Matrix Cheb3D(double t, int N, double Ta, double Tb, Matrix& Cx, Matrix& Cy, Matrix& Cz);


#endif // CHEB3D_HPP