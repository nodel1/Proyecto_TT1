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

void IERS(double Mjd_UTC, char interp, double& x_pole, double& y_pole, double& UT1_UTC, double& LOD, double& dpsi, double& deps, double& dx_pole, double& dy_pole, double& TAI_UTC) {
    double x_pole;
    double y_pole;
    double UT1_UTC;
    double LOD;
    double dpsi;
    double deps;
    double dx_pole;
    double dy_pole;
    double TAI_UTC;
    Matrix& eop = eopdata;

    if (interp == 'l') {
        double mjd = floor(Mjd_UTC);
        int i = -1;

        for (int col = 0; col < eop.n_column; col++) {
            if (mjd == eop(3, col)) {
                i = col;
                break;
            }
        }
        if (i == -1) {
            cout << "IERS: MJD not found in eop data.\n";
            exit(EXIT_FAILURE);
        }

        Matrix preeop = zeros(13, 1);
        Matrix nexteop = zeros(13, 1);
        for (int k = 0; k < 13; k++) {
            preeop(k, 0) = eop(k, i);
            nexteop(k, 0) = eop(k, i + 1);
        }

        double mfme = 1440.0 * (Mjd_UTC - mjd);
        double fixf = mfme / 1440.0;

        // Setting of IERS Earth rotation parameters
        // (UT1-UTC [s], TAI-UTC [s], x ["], y ["])
        x_pole = preeop(4, 0) + (nexteop(4, 0) - preeop(4, 0)) * fixf;
        y_pole = preeop(5, 0) + (nexteop(5, 0) - preeop(5, 0)) * fixf;
        UT1_UTC = preeop(6, 0) + (nexteop(6, 0) - preeop(6, 0)) * fixf;
        LOD = preeop(7, 0) + (nexteop(7, 0) - preeop(7, 0)) * fixf;
        dpsi = preeop(8, 0) + (nexteop(8, 0) - preeop(8, 0)) * fixf;
        deps = preeop(9, 0) + (nexteop(9, 0) - preeop(9, 0)) * fixf;
        dx_pole = preeop(10, 0) + (nexteop(10, 0) - preeop(10, 0)) * fixf;
        dy_pole = preeop(11, 0) + (nexteop(11, 0) - preeop(11, 0)) * fixf;
        TAI_UTC = preeop(12, 0);

        x_pole /= Const::Arcs;  // Pole coordinate [rad]
        y_pole /= Const::Arcs;  // Pole coordinate [rad]
        dpsi /= Const::Arcs;
        deps /= Const::Arcs;
        dx_pole /= Const::Arcs; // Pole coordinate [rad]
        dy_pole /= Const::Arcs; // Pole coordinate [rad]
    } else if (interp == 'n') {
        double mjd = floor(Mjd_UTC);
        int i = -1;

        for (int col = 0; col < eop.n_column; col++) {
            if (mjd == eop(3, col)) {
                i = col;
                break;
            }
        }
        if (i == -1) {
            cout << "IERS: MJD not found in eop data.\n";
            exit(EXIT_FAILURE);
        }

        // Setting of IERS Earth rotation parameters
        // (UT1-UTC [s], TAI-UTC [s], x ["], y ["])
        x_pole = eop(4, i) / Const::Arcs;  // Pole coordinate [rad]
        y_pole = eop(5, i) / Const::Arcs;  // Pole coordinate [rad]
        UT1_UTC = eop(6, i);               // UT1-UTC time difference [s]
        LOD = eop(7, i);                   // Length of day [s]
        dpsi = eop(8, i) / Const::Arcs;
        deps = eop(9, i) / Const::Arcs;
        dx_pole = eop(10, i) / Const::Arcs; // Pole coordinate [rad]
        dy_pole = eop(11, i) / Const::Arcs; // Pole coordinate [rad]
        TAI_UTC = eop(12, i);               // TAI-UTC time difference [s]
    } else {
        cout << "IERS: Invalid interpolation method '" << interp << "'\n";
        exit(EXIT_FAILURE);
    }
}

void IERS(double Mjd_UTC, double& x_pole, double& y_pole, double& UT1_UTC, double& LOD, double& dpsi, double& deps, double& dx_pole, double& dy_pole, double& TAI_UTC) {
    IERS(Mjd_UTC, 'n', x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC);
}