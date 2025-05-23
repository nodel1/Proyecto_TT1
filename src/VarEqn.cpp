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
/** @file VarEqn.cpp
 *  @brief Computes variational equations for orbit propagation
 *
 *  @author Noel Del Rio
 *  @bug No known bugs
 */

#include "../include/VarEqn.hpp"
#include "../include/IERS.hpp"
#include "../include/timediff.hpp"
#include "../include/PrecMatrix.hpp"
#include "../include/NutMatrix.hpp"
#include "../include/PoleMatrix.hpp"
#include "../include/GHAMatrix.hpp"
#include "../include/AccelHarmonic.hpp"
#include "../include/G_AccelHarmonic.hpp"
#include "../include/global.hpp"

Matrix VarEqn(double x, Matrix yPhi) {
    double x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC;
    double Mjd_UTC = AuxParam.Mjd_UTC + x / 86400.0;

    IERS(Mjd_UTC, 'l', x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC);

    double UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC;
    timediff(UT1_UTC, TAI_UTC, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);

    double Mjd_UT1 = Mjd_UTC + UT1_UTC / 86400.0;
    double Mjd_TT = Mjd_UTC + TT_UTC / 86400.0;

    Matrix P = PrecMatrix(Const::MJD_J2000, Mjd_TT);
    Matrix N = NutMatrix(Mjd_TT);
    Matrix T = N * P;

    Matrix GHA = GHAMatrix(Mjd_UT1);
	Matrix temp = GHA * T;
	Matrix E = PoleMatrix(x_pole, y_pole) * temp;

    Matrix r(3, 1), v(3, 1);
    for (int i = 1; i <= 3; i++) {
        r(i, 1) = yPhi(i, 1);
        v(i, 1) = yPhi(i + 3, 1);
    }

    Matrix Phi(6, 6);
    for (int j = 1; j <= 6; j++)
        for (int i = 1; i <= 6; i++)
            Phi(i, j) = yPhi(6 * (j - 1) + i + 6, 1);

    Matrix a = AccelHarmonic(r, E, AuxParam.n, AuxParam.m);
    Matrix G = G_AccelHarmonic(r, E, AuxParam.n, AuxParam.m);

    Matrix dfdy = zeros(6, 6);
    for (int i = 1; i <= 3; ++i) {
        for (int j = 1; j <= 3; ++j) {
            dfdy(i + 3, j) = G(i, j);
            dfdy(i, j + 3) = (i == j) ? 1.0 : 0.0;
        }
    }

    Matrix Phip = dfdy * Phi;
    Matrix yPhip(42, 1);

    for (int i = 1; i <= 3; i++) {
        yPhip(i, 1) = v(i, 1);
        yPhip(i + 3, 1) = a(i, 1);
    }

    for (int j = 1; j <= 6; j++)
        for (int i = 1; i <= 6; i++)
            yPhip(6 * (j - 1) + i + 6, 1) = Phip(i, j);

    return yPhip;
}
