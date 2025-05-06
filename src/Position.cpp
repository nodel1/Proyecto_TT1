// $Source$
//--------------------------------------------------
// Position
//--------------------------------------------------
// Proyecto_TT1: Position Vector Calculation Implementation
//
// Copyright (c) 2025, Meysam Mahooti
//
// Created: 2025/04/23
//
/** @file Position.cpp
 *  @brief Implementation of position vector calculation from geodetic coordinates
 *
 *  @author Noel Del Rio Gonzalez
 *  @bug No known bugs
 */

#include "..\include\Position.hpp"
#include "..\include\SAT_Const.hpp"

using namespace std;

Matrix Position(double lon,double  lat, double h){
	
	double R_equ = Const::R_Earth;
	double f     = Const::f_Earth;


	double e2     = f*(2.0-f);   
	double CosLat = cos(lat);    
	double SinLat = sin(lat);

	double N = R_equ / sqrt(1.0-e2*SinLat*SinLat);
	
    Matrix r(3, 1);
    r(1, 1) = (N + h) * CosLat * cos(lon);
    r(2, 1) = (N + h) * CosLat * sin(lon);
    r(3, 1) = ((1.0 - e2) * N + h) * SinLat;
	
	return r;
}