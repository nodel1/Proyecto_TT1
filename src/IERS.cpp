// $Source$
//--------------------------------------------------
// IERS
//--------------------------------------------------
// Proyecto_TT1: IERS Time and Polar Motion Data Management Implementation
//
// Copyright (c) 2025, Meysam Mahooti
//
// Created: 2025/04/23
//
/** @file IERS.cpp
 *  @brief Implementation of IERS time and polar motion data management
 *
 *  @author Noel Del Rio Gonzalez
 *  @bug No known bugs
 */

#include "../include/IERS.hpp"

using namespace std;

void IERS(double Mjd_UTC, char interp, double& x_pole, double& y_pole, double& UT1_UTC, double& LOD, double& dpsi, double& deps, double& dx_pole, double& dy_pole, double& TAI_UTC) {
    Matrix& eop = eopdata; // Variable global definida en global.hpp

    if (eop.n_row < 13) {
        cout << "IERS: eop matrix must have at least 13 rows.\n";
        exit(EXIT_FAILURE);
    }

    if (interp == 'l') {
        // Interpolación lineal
        double mjd = floor(Mjd_UTC);
        int i = -1;

        // Buscar el índice donde eop(4, col) == mjd
        for (int col = 1; col <= eop.n_column; col++) {
            if (fabs(mjd - eop(4, col)) < 1e-10) { // Tolerancia para comparación de doubles
                i = col;
                break;
            }
        }
        if (i == -1 || i >= eop.n_column) {
            cout << "IERS: MJD " << mjd << " not found in eop data or insufficient data for interpolation.\n";
            exit(EXIT_FAILURE);
        }

        // Extraer columnas preeop y nexteop
        Matrix preeop = extract_column(eop, i);
        Matrix nexteop = extract_column(eop, i + 1);

        double mfme = 1440.0 * (Mjd_UTC - mjd);
        double fixf = mfme / 1440.0;

        // Interpolación lineal de los parámetros
        x_pole  = preeop(5, 1) + (nexteop(5, 1) - preeop(5, 1)) * fixf;
        y_pole  = preeop(6, 1) + (nexteop(6, 1) - preeop(6, 1)) * fixf;
        UT1_UTC = preeop(7, 1) + (nexteop(7, 1) - preeop(7, 1)) * fixf;
        LOD     = preeop(8, 1) + (nexteop(8, 1) - preeop(8, 1)) * fixf;
        dpsi    = preeop(9, 1) + (nexteop(9, 1) - preeop(9, 1)) * fixf;
        deps    = preeop(10, 1) + (nexteop(10, 1) - preeop(10, 1)) * fixf;
        dx_pole = preeop(11, 1) + (nexteop(11, 1) - preeop(11, 1)) * fixf;
        dy_pole = preeop(12, 1) + (nexteop(12, 1) - preeop(12, 1)) * fixf;
        TAI_UTC = preeop(13, 1); // TAI_UTC no se interpola

        // Convertir de arcosegundos a radianes
        x_pole  /= Const::Arcs;
        y_pole  /= Const::Arcs;
        dpsi    /= Const::Arcs;
        deps    /= Const::Arcs;
        dx_pole /= Const::Arcs;
        dy_pole /= Const::Arcs;
    } else if (interp == 'n') {
        // Sin interpolación
        double mjd = floor(Mjd_UTC);
        int i = -1;

        // Buscar el índice donde eop(4, col) == mjd
        for (int col = 1; col <= eop.n_column; col++) {
            if (fabs(mjd - eop(4, col)) < 1e-10) {
                i = col;
                break;
            }
        }
        if (i == -1) {
            cout << "IERS: MJD " << mjd << " not found in eop data.\n";
            exit(EXIT_FAILURE);
        }

        // Asignar parámetros directamente
        x_pole  = eop(5, i) / Const::Arcs;
        y_pole  = eop(6, i) / Const::Arcs;
        UT1_UTC = eop(7, i);
        LOD     = eop(8, i);
        dpsi    = eop(9, i) / Const::Arcs;
        deps    = eop(10, i) / Const::Arcs;
        dx_pole = eop(11, i) / Const::Arcs;
        dy_pole = eop(12, i) / Const::Arcs;
        TAI_UTC = eop(13, i);
    } else {
        cout << "IERS: Invalid interpolation method '" << interp << "'\n";
        exit(EXIT_FAILURE);
    }
}

void IERS(double Mjd_UTC, double& x_pole, double& y_pole, double& UT1_UTC, double& LOD, double& dpsi, double& deps, double& dx_pole, double& dy_pole, double& TAI_UTC) {
    IERS(Mjd_UTC, 'n', x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC);
}