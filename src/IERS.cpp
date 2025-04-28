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

#include "..\include\IERS.hpp"

using namespace std;

// Helper function to find the first index where eop(3,j) == mjd
int find_index(const Matrix& eop, double mjd) {
    for (int j = 0; j < eop.n_column; j++) {
        if (eop(3, j) == mjd) {
            return j;
        }
    }
    cout << "Error: MJD " << mjd << " not found in eop data" << endl;
    exit(EXIT_FAILURE);
    return -1; // Unreachable, but included for completeness
}

void IERS(const Matrix& eop, double Mjd_UTC, char interp, double& x_pole, double& y_pole, double& UT1_UTC, double& LOD, double& dpsi, double& deps, double& dx_pole, double& dy_pole, double& TAI_UTC) {
    double mjd = floor(Mjd_UTC);
    int i = find_index(eop, mjd);

    if (interp == 'l') {
        // Linear interpolation
        Matrix preeop(13, 1);
        Matrix nexteop(13, 1);
        for (int k = 0; k < 13; k++) {
            preeop(k, 0) = eop(k, i);
            nexteop(k, 0) = eop(k, i + 1);
        }

        double mfme = 1440.0 * (Mjd_UTC - mjd);
        double fixf = mfme / 1440.0;

        // Setting of IERS Earth rotation parameters
        x_pole = preeop(4, 0) + (nexteop(4, 0) - preeop(4, 0)) * fixf;
        y_pole = preeop(5, 0) + (nexteop(5, 0) - preeop(5, 0)) * fixf;
        UT1_UTC = preeop(6, 0) + (nexteop(6, 0) - preeop(6, 0)) * fixf;
        LOD = preeop(7, 0) + (nexteop(7, 0) - preeop(7, 0)) * fixf;
        dpsi = preeop(8, 0) + (nexteop(8, 0) - preeop(8, 0)) * fixf;
        deps = preeop(9, 0) + (nexteop(9, 0) - preeop(9, 0)) * fixf;
        dx_pole = preeop(10, 0) + (nexteop(10, 0) - preeop(10, 0)) * fixf;
        dy_pole = preeop(11, 0) + (nexteop(11, 0) - preeop(11, 0)) * fixf;
        TAI_UTC = preeop(12, 0);

        // Convert to radians
        x_pole /= Const::Arcs;
        y_pole /= Const::Arcs;
        dpsi /= Const::Arcs;
        deps /= Const::Arcs;
        dx_pole /= Const::Arcs;
        dy_pole /= Const::Arcs;
    } else if (interp == 'n') {
        // No interpolation
        x_pole = eop(4, i) / Const::Arcs;
        y_pole = eop(5, i) / Const::Arcs;
        UT1_UTC = eop(6, i);
        LOD = eop(7, i);
        dpsi = eop(8, i) / Const::Arcs;
        deps = eop(9, i) / Const::Arcs;
        dx_pole = eop(10, i) / Const::Arcs;
        dy_pole = eop(11, i) / Const::Arcs;
        TAI_UTC = eop(12, i);
    } else {
        cout << "Error: Invalid interpolation method '" << interp << "'" << endl;
        exit(EXIT_FAILURE);
    }
}

void IERS(const Matrix& eop, double Mjd_UTC, double& x_pole, double& y_pole, double& UT1_UTC, double& LOD, double& dpsi, double& deps, double& dx_pole, double& dy_pole, double& TAI_UTC) {
    IERS(eop, Mjd_UTC, 'n', x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC);
}