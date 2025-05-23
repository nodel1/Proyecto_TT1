// $Source$
//--------------------------------------------------
// MeasUpdate
//--------------------------------------------------
// Proyecto_TT1: Kalman Filter Measurement Update Implementation
//
// Copyright (c) 2025, Meysam Mahooti
//
// Created: 2025/05/21
//
/** @file MeasUpdate.cpp
 *  @brief Implementation of Kalman filter measurement update step for scalar measurements
 *
 *  @author Noel Del Rio Gonzalez
 *  @bug No known bugs
 */

#include "..\include\MeasUpdate.hpp"
#include <cmath>
#include <iostream>

using namespace std;

void MeasUpdate(Matrix& x, double z, double g, double s, 
                 Matrix& G, Matrix& P, int n, Matrix& K)
{
    // Kalman gain calculation (scalar measurement version)
	    cout << "llegas aqui1" << endl;
		
    Matrix G_transpose = transpose(G);
	
		    cout << "llegas aqui2" << endl;
	
    double S = s*s + (G*P*G_transpose)(1,1);  // Innovation covariance (scalar)
	
		    cout << "llegas aqui3" << endl;
			
			
    K = (P*G_transpose)*(1.0/S);  // Kalman gain (n x 1)


	    cout << "llegas aqui4" << endl;
		
		
    // State update
    double innovation = z - g;
	
		    cout << "llegas aqui5" << endl;
			
			
    x = x + K*innovation;
	
		    cout << "llegas aqui6" << endl;
			

    // Covariance update (Joseph form for better numerical stability)
    Matrix I = eye(n);
	
		    cout << "llegas aqui7" << endl;
			
			
    Matrix I_KG = I - K*G;
	
		    cout << "llegas aqui8" << endl;
			
			
    P = I_KG*P*(transpose(I_KG)) + K*(s*s)*transpose(K);
	
		    cout << "llegas aqui final" << endl;
			
}