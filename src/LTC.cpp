// $Source$
//--------------------------------------------------
// LTC
//--------------------------------------------------
// Proyecto_TT1: Local Tangent Coordinates Transformation
//
// Copyright (c) 2025, Meysam Mahooti
//
// Created: 2025/05/07
//
/** @file LTC.cpp
 *  @brief Transformation from Greenwich meridian system to local tangent coordinates
 *
 *  @author Noel Del Rio Gonzalez
 *  @bug No known bugs
 */

#include "..\include\LTC.hpp"
#include "..\include\R_y.hpp"
#include "..\include\R_z.hpp"

using namespace std;

Matrix LTC(double lon, double lat) {
    // Calcular la matriz de rotación: M = R_y(-lat) * R_z(lon)
    Matrix M = R_y(-lat) * R_z(lon);

    // Permutación cíclica de las filas: 1->2, 2->3, 3->1
    Matrix temp = zeros(3, 3); // Matriz temporal para almacenar el resultado
    for (int j = 1; j <= 3; j++) {
        temp(1, j) = M(2, j); // Fila 2 de M -> Fila 1 de temp
        temp(2, j) = M(3, j); // Fila 3 de M -> Fila 2 de temp
        temp(3, j) = M(1, j); // Fila 1 de M -> Fila 3 de temp
    }

    return temp;
}