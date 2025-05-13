#include "..\include\matrix.hpp"
#include <cstdio>
#include <cmath>

#include "..\include\SAT_Const.hpp"
#include "..\include\global.hpp"

#include "..\include\AccelPointMass.hpp"
#include "..\include\Cheb3D.hpp"
#include "..\include\EccAnom.hpp"
#include "..\include\Frac.hpp"
#include "..\include\MeanObliquity.hpp"
#include "..\include\Mjday.hpp"
#include "..\include\Mjday_TDB.hpp"
#include "..\include\Position.hpp"
#include "..\include\R_x.hpp"
#include "..\include\R_y.hpp"
#include "..\include\R_z.hpp"
#include "..\include\sign_.hpp"
#include "..\include\timediff.hpp"
#include "..\include\AzElPa.hpp"
#include "..\include\IERS.hpp"
#include "..\include\Legendre.hpp"
#include "..\include\NutAngles.hpp"
#include "..\include\TimeUpdate.hpp"
#include "..\include\LTC.hpp"
#include "..\include\NutMatrix.hpp"
#include "..\include\PoleMatrix.hpp"
#include "..\include\PrecMatrix.hpp"
#include "..\include\gmst.hpp"
#include "..\include\AccelHarmonic.hpp"
#include "..\include\EqnEquinox.hpp"
#include "..\include\JPL_Eph_DE430.hpp"


using namespace std;




int tests_run = 0;

#define FAIL() printf("\nfailure in %s() line %d\n", __func__, __LINE__)
#define _assert(test) do { if (!(test)) { FAIL(); return 1; } } while(0)
#define _verify(test) do { int r=test(); tests_run++; if(r) return r; } while(0)

int m_equals(Matrix A, Matrix B, double p) {
    if (A.n_row != B.n_row || A.n_column != B.n_column)
        return 0;
    else
        for(int i = 1; i <= A.n_row; i++)
            for(int j = 1; j <= A.n_column; j++)
                if(fabs(A(i,j)-B(i,j)) > p) {
                    printf("%2.20lf %2.20lf\n",A(i,j),B(i,j));
                    return 0;
                }

    return 1;
}

// Test para Matrix(row, column)
int m_constructor_01() {
    Matrix A(2, 3);
    A(1,1) = 1; A(1,2) = 2; A(1,3) = 3;
    A(2,1) = 4; A(2,2) = 5; A(2,3) = 6;

    Matrix expected(2, 3);
    expected(1,1) = 1; expected(1,2) = 2; expected(1,3) = 3;
    expected(2,1) = 4; expected(2,2) = 5; expected(2,3) = 6;

    _assert(m_equals(A, expected, 1e-10));
    return 0;
}

// Test para operator () (row, column)
int m_access_01() {
    Matrix A(2, 2);
    A(1,1) = 1; A(1,2) = 2;
    A(2,1) = 3; A(2,2) = 4;

    _assert(fabs(A(1,1) - 1) < 1e-10);
    _assert(fabs(A(1,2) - 2) < 1e-10);
    _assert(fabs(A(2,1) - 3) < 1e-10);
    _assert(fabs(A(2,2) - 4) < 1e-10);
    return 0;
}

// Test para operator + (Matrix)
int m_sum_01() {
    int f = 3;
    int c = 4;

    Matrix A(f, c);
    A(1,1) = 0; A(1,2) =  2; A(1,3) = 8; A(1,4) = 0;
    A(2,1) = 1; A(2,2) = -1; A(2,3) = 0; A(2,4) = 0;
    A(3,1) = 0; A(3,2) =  1; A(3,3) = 0; A(3,4) = 5;

    Matrix B(f, c);
    B(1,1) = 2; B(1,2) =  0; B(1,3) = 0; B(1,4) = 0;
    B(2,1) = 7; B(2,2) = -2; B(2,3) = 1; B(2,4) = 0;
    B(3,1) = 0; B(3,2) = -3; B(3,3) = 0; B(3,4) = 2;

    Matrix C(f, c);
    C(1,1) = 2; C(1,2) =  2; C(1,3) = 8; C(1,4) = 0;
    C(2,1) = 8; C(2,2) = -3; C(2,3) = 1; C(2,4) = 0;
    C(3,1) = 0; C(3,2) = -2; C(3,3) = 0; C(3,4) = 7;

    Matrix R = A + B;

    _assert(m_equals(C, R, 1e-10));
    return 0;
}

// Test para operator - (Matrix)
int m_sub_01() {
    int f = 3;
    int c = 4;

    Matrix A(f, c);
    A(1,1) = 0; A(1,2) = 2; A(1,3) = 8; A(1,4) = 0;
    A(2,1) = 1; A(2,2) = -1; A(2,3) = 0; A(2,4) = 0;
    A(3,1) = 0; A(3,2) = 1; A(3,3) = 0; A(3,4) = 5;

    Matrix B(f, c);
    B(1,1) = 2; B(1,2) = 0; B(1,3) = 0; B(1,4) = 0;
    B(2,1) = 7; B(2,2) = -2; B(2,3) = 1; B(2,4) = 0;
    B(3,1) = 0; B(3,2) = -3; B(3,3) = 0; B(3,4) = 2;

    Matrix C(f, c);
    C(1,1) = -2; C(1,2) = 2; C(1,3) = 8; C(1,4) = 0;
    C(2,1) = -6; C(2,2) = 1; C(2,3) = -1; C(2,4) = 0;
    C(3,1) = 0; C(3,2) = 4; C(3,3) = 0; C(3,4) = 3;

    Matrix R = A - B;

    _assert(m_equals(C, R, 1e-10));
    return 0;
}

// Test para operator * (Matrix)
int m_mult_01() {
    Matrix A(2, 3);
    A(1,1) = 1; A(1,2) = 2; A(1,3) = 3;
    A(2,1) = 4; A(2,2) = 5; A(2,3) = 6;

    Matrix B(3, 2);
    B(1,1) = 7; B(1,2) = 8;
    B(2,1) = 9; B(2,2) = 10;
    B(3,1) = 11; B(3,2) = 12;

    Matrix C(2, 2);
    C(1,1) = 58; C(1,2) = 64;
    C(2,1) = 139; C(2,2) = 154;

    Matrix R = A * B;

    _assert(m_equals(C, R, 1e-10));
    return 0;
}

// Test para operator / (Matrix)
int m_div_01() {
    Matrix A(2, 2);
    A(1,1) = 1; A(1,2) = 2;
    A(2,1) = 3; A(2,2) = 4;

    Matrix B(2, 2);
    B(1,1) = 4; B(1,2) = 7;
    B(2,1) = 2; B(2,2) = 6;

    Matrix expected(2, 2);
    expected(1,1) = 0.2; expected(1,2) = 0.1;
    expected(2,1) = 1.0; expected(2,2) = -0.5;

    Matrix R = A / B;

    _assert(m_equals(expected, R, 1e-10));
    return 0;
}

// Test para operator = (Matrix)
int m_assign_01() {
    Matrix A(2, 2);
    A(1,1) = 1; A(1,2) = 2;
    A(2,1) = 3; A(2,2) = 4;

    Matrix B(2, 2);
    B = A;

    Matrix expected(2, 2);
    expected(1,1) = 1; expected(1,2) = 2;
    expected(2,1) = 3; expected(2,2) = 4;

    _assert(m_equals(B, expected, 1e-10));
    return 0;
}

// Test para zeros(row, column)
int m_zeros_01() {
    int f = 3;
    int c = 4;

    Matrix A(f, c);
    A(1,1) = 0; A(1,2) = 0; A(1,3) = 0; A(1,4) = 0;
    A(2,1) = 0; A(2,2) = 0; A(2,3) = 0; A(2,4) = 0;
    A(3,1) = 0; A(3,2) = 0; A(3,3) = 0; A(3,4) = 0;

    Matrix B = zeros(3, 4);

    _assert(m_equals(A, B, 1e-10));
    return 0;
}

// Test para eye(n)
int m_eye_01() {
    Matrix expected(3, 3);
    expected(1,1) = 1; expected(1,2) = 0; expected(1,3) = 0;
    expected(2,1) = 0; expected(2,2) = 1; expected(2,3) = 0;
    expected(3,1) = 0; expected(3,2) = 0; expected(3,3) = 1;

    Matrix R = eye(3);

    _assert(m_equals(expected, R, 1e-10));
    return 0;
}

// Test para transpose(Matrix)
int m_transpose_01() {
    Matrix A(2, 3);
    A(1,1) = 1; A(1,2) = 2; A(1,3) = 3;
    A(2,1) = 4; A(2,2) = 5; A(2,3) = 6;

    Matrix expected(3, 2);
    expected(1,1) = 1; expected(1,2) = 4;
    expected(2,1) = 2; expected(2,2) = 5;
    expected(3,1) = 3; expected(3,2) = 6;

    Matrix R = transpose(A);

    _assert(m_equals(expected, R, 1e-10));
    return 0;
}

// Test para inv(Matrix)
int m_inv_01() {
    Matrix A(2, 2);
    A(1,1) = 4; A(1,2) = 7;
    A(2,1) = 2; A(2,2) = 6;

    Matrix B(2, 2);
    B(1,1) = 0.6; B(1,2) = -0.7;
    B(2,1) = -0.2; B(2,2) = 0.4;

    Matrix R = inv(A);

    _assert(m_equals(B, R, 1e-10));
    return 0;
}

// Test para operator + (double)
int m_scalar_add_01() {
    Matrix A(2, 2);
    A(1,1) = 1; A(1,2) = 2;
    A(2,1) = 3; A(2,2) = 4;

    Matrix expected(2, 2);
    expected(1,1) = 3; expected(1,2) = 4;
    expected(2,1) = 5; expected(2,2) = 6;

    Matrix R = A + 2.0;

    _assert(m_equals(expected, R, 1e-10));
    return 0;
}

// Test para operator - (double)
int m_scalar_sub_01() {
    Matrix A(2, 2);
    A(1,1) = 3; A(1,2) = 4;
    A(2,1) = 5; A(2,2) = 6;

    Matrix expected(2, 2);
    expected(1,1) = 1; expected(1,2) = 2;
    expected(2,1) = 3; expected(2,2) = 4;

    Matrix R = A - 2.0;

    _assert(m_equals(expected, R, 1e-10));
    return 0;
}

// Test para operator * (double)
int m_scalar_mult_01() {
    Matrix A(2, 2);
    A(1,1) = 1; A(1,2) = 2;
    A(2,1) = 3; A(2,2) = 4;

    Matrix expected(2, 2);
    expected(1,1) = 2; expected(1,2) = 4;
    expected(2,1) = 6; expected(2,2) = 8;

    Matrix R = A * 2.0;

    _assert(m_equals(expected, R, 1e-10));
    return 0;
}

// Test para operator / (double)
int m_scalar_div_01() {
    Matrix A(2, 2);
    A(1,1) = 2; A(1,2) = 4;
    A(2,1) = 6; A(2,2) = 8;

    Matrix expected(2, 2);
    expected(1,1) = 1; expected(1,2) = 2;
    expected(2,1) = 3; expected(2,2) = 4;

    Matrix R = A / 2.0;

    _assert(m_equals(expected, R, 1e-10));
    return 0;
}

// Test para Matrix(n)
int m_constructor_n_01() {
    Matrix A(2);
    A(1,1) = 1; A(1,2) = 2;
    A(2,1) = 3; A(2,2) = 4;

    Matrix expected(2, 2);
    expected(1,1) = 1; expected(1,2) = 2;
    expected(2,1) = 3; expected(2,2) = 4;

    _assert(m_equals(A, expected, 1e-10));
    return 0;
}

// Test para operator () (n)
int m_vector_access_01() {
    Matrix A(1, 3);
    A(1,1) = 1; A(1,2) = 2; A(1,3) = 3;

    _assert(fabs(A(1) - 1) < 1e-10);
    _assert(fabs(A(2) - 2) < 1e-10);
    _assert(fabs(A(3) - 3) < 1e-10);
    return 0;
}

// Test para zeros(n)
int m_zeros_n_01() {
    Matrix expected(2, 2);
    expected(1,1) = 0; expected(1,2) = 0;
    expected(2,1) = 0; expected(2,2) = 0;

    Matrix R = zeros(2);

    _assert(m_equals(expected, R, 1e-10));
    return 0;
}

// Test para norm
int m_norm_01() {
    Matrix A(2, 2);
    A(1,1) = 1; A(1,2) = 2;
    A(2,1) = 3; A(2,2) = 4;

    double expected = sqrt(1*1 + 2*2 + 3*3 + 4*4); // sqrt(30)
    double result = norm(A);

    _assert(fabs(expected - result) < 1e-10);
    return 0;
}

// Test para dot
int m_dot_01() {
    Matrix A(1, 3);
    A(1,1) = 1; A(1,2) = 2; A(1,3) = 3;

    Matrix B(1, 3);
    B(1,1) = 4; B(1,2) = 5; B(1,3) = 6;

    double expected = 1*4 + 2*5 + 3*6; // 32
    double result = dot(A, B);

    _assert(fabs(expected - result) < 1e-10);
    return 0;
}

// Test para cross
int m_cross_01() {
    Matrix A(3, 1);
    A(1,1) = 1; A(2,1) = 2; A(3,1) = 3;

    Matrix B(3, 1);
    B(1,1) = 4; B(2,1) = 5; B(3,1) = 6;

    Matrix expected(3, 1);
    expected(1,1) = 2*6 - 3*5; // -3
    expected(2,1) = 3*4 - 1*6; // 6
    expected(3,1) = 1*5 - 2*4; // -3

    Matrix R = cross(A, B);

    _assert(m_equals(expected, R, 1e-10));
    return 0;
}

// Test para extract_vector
int m_extract_vector_01() {
    Matrix A(2, 2);
    A(1,1) = 1; A(1,2) = 2;
    A(2,1) = 3; A(2,2) = 4;

    Matrix expected(1, 4);
    expected(1,1) = 1; expected(1,2) = 2; expected(1,3) = 3; expected(1,4) = 4;

    Matrix R = extract_vector(A);

    _assert(m_equals(expected, R, 1e-10));
    return 0;
}

// Test para union_vector
int m_union_vector_01() {
    Matrix V(1, 4);
    V(1,1) = 1; V(1,2) = 2; V(1,3) = 3; V(1,4) = 4;

    Matrix expected(2, 2);
    expected(1,1) = 1; expected(1,2) = 2;
    expected(2,1) = 3; expected(2,2) = 4;

    Matrix R = union_vector(V, 2, 2);

    _assert(m_equals(expected, R, 1e-10));
    return 0;
}

// Test para extract_row
int m_extract_row_01() {
    Matrix A(2, 3);
    A(1,1) = 1; A(1,2) = 2; A(1,3) = 3;
    A(2,1) = 4; A(2,2) = 5; A(2,3) = 6;

    Matrix expected(1, 3);
    expected(1,1) = 4; expected(1,2) = 5; expected(1,3) = 6;

    Matrix R = extract_row(A, 2);

    _assert(m_equals(expected, R, 1e-10));
    return 0;
}

// Test para extract_column
int m_extract_column_01() {
    Matrix A(2, 3);
    A(1,1) = 1; A(1,2) = 2; A(1,3) = 3;
    A(2,1) = 4; A(2,2) = 5; A(2,3) = 6;

    Matrix expected(2, 1);
    expected(1,1) = 2; expected(2,1) = 5;

    Matrix R = extract_column(A, 2);

    _assert(m_equals(expected, R, 1e-10));
    return 0;
}

// Test para assign_row
int m_assign_row_01() {
    Matrix A(2, 3);
    A(1,1) = 1; A(1,2) = 2; A(1,3) = 3;
    A(2,1) = 4; A(2,2) = 5; A(2,3) = 6;

    Matrix V(1, 3);
    V(1,1) = 7; V(1,2) = 8; V(1,3) = 9;

    assign_row(A, 1, V);

    Matrix expected(2, 3);
    expected(1,1) = 7; expected(1,2) = 8; expected(1,3) = 9;
    expected(2,1) = 4; expected(2,2) = 5; expected(2,3) = 6;

    _assert(m_equals(A, expected, 1e-10));
    return 0;
}

// Test para assign_column
int m_assign_column_01() {
    Matrix A(2, 3);
    A(1,1) = 1; A(1,2) = 2; A(1,3) = 3;
    A(2,1) = 4; A(2,2) = 5; A(2,3) = 6;

    Matrix V(2, 1);
    V(1,1) = 7; V(2,1) = 8;

    assign_column(A, 2, V);

    Matrix expected(2, 3);
    expected(1,1) = 1; expected(1,2) = 7; expected(1,3) = 3;
    expected(2,1) = 4; expected(2,2) = 8; expected(2,3) = 6;

    _assert(m_equals(A, expected, 1e-10));
    return 0;
}


// Test AccelPointMass
int m_accel_point_mass_01() {
    Matrix r(1, 3);
    r(1,1) = Const::R_Earth + 700e3; // 7078136.3 m
    r(1,2) = 0;
    r(1,3) = 0;

    Matrix s(1, 3);
    s(1,1) = 384400e3; // 384400000 m
    s(1,2) = 0;
    s(1,3) = 0;

    double GM = Const::GM_Moon; // 4.9028e12 m^3/s^2

    Matrix expected(1, 3);
    expected(1,1) = 1.25651832363664e-06;
    expected(1,2) = 0;
    expected(1,3) = 0;

    Matrix a = AccelPointMass(r, s, GM);

    _assert(m_equals(a, expected, 1e-9));     //resultados sacados de hacerpruebas en matlab

    double expected_magnitude = 1.25651832363664e-06;
    double result_magnitude = norm(a);
    _assert(fabs(expected_magnitude - result_magnitude) < 1e-9);

    return 0;
}

// Test Cheb3D
int m_cheb3d_01() {
    int N = 3;
    double Ta = 0.0;
    double Tb = 1.0;
    double t = 0.5;
    Matrix Cx(1, 3);
    Cx(1, 1) = 1.0; Cx(1, 2) = 0.5; Cx(1, 3) = 0.3;
    Matrix Cy(1, 3);
    Cy(1, 1) = 0.0; Cy(1, 2) = 0.1; Cy(1, 3) = 0.2;
    Matrix Cz(1, 3);
    Cz(1, 1) = 0.5; Cz(1, 2) = 0.1; Cz(1, 3) = 0.2;
    Matrix R = Cheb3D(t, N, Ta, Tb, Cx, Cy, Cz);
    Matrix expected(1, 3);
    expected(1, 1) = 0.7; expected(1, 2) = -0.2; expected(1, 3) = 0.3;
    _assert(m_equals(expected, R, 1e-8));
    return 0;
}


// Test EccAnom
int m_ecc_anom_01() {     //pruebo con diferentes valores
    // (e=0)
    _assert(fabs(EccAnom(0.7854, 0.0) - 0.7854) < 1e-4);

    //(e=0.1)
    _assert(fabs(EccAnom(1.5708, 0.1) - 1.6703) < 1e-4);

    // (e=0.5)
    _assert(fabs(EccAnom(3.1416, 0.5) - 3.1416) < 1e-4);

    // (e=0.9)
    _assert(fabs(EccAnom(4.7124, 0.9) - 4.0198) < 1e-4);

    return 0;      //VALORES DE PROBAR EN MATLAB
}

int m_frac_01() {
    double x = 3.75;
    double expected = 0.75;
    double result = Frac(x);


    _assert(fabs(result - expected) < 1e-10);
    return 0;
}

// Test MeanObliquity
int m_mean_obliquity_01() {

    double Mjd_TT1 = Const::MJD_J2000;
    double expected1 = 0.409092804222;
    double result1 = MeanObliquity(Mjd_TT1);
    _assert(fabs(result1 - expected1) < 1e-10);


    double Mjd_TT2 = Const::MJD_J2000 + 36525;
    double expected2 = 0.408865844627;
    double result2 = MeanObliquity(Mjd_TT2);
    _assert(fabs(result2 - expected2) < 1e-10);

    return 0;
}

//test MJDay
int m_mjday_01() {
    double result = Mjday(2025, 1, 1, 0, 0, 0.0);
    double expected = 60676.0;
    _assert(fabs(result - expected) < 1e-9);
    return 0;
}

//test Mjday_TDB
int m_mjday_tdb_01() {
    double result = Mjday_TDB(60355.0);
    double expected = 60355.000000012;
    _assert(fabs(result - expected) < 1e-9);
    return 0;
}

//test position
int m_position_01() {

    Matrix r = Position(0.0, 0.0, 0.0);


    Matrix expected(3, 1);
    expected(1, 1) = 6378136.3;
    expected(2, 1) = 0.0;
    expected(3, 1) = 0.0;

    _assert(m_equals(r, expected, 1e-6));  // Precisión tolerada
    return 0;
}

//test R_x
int m_r_x_01() {
    Matrix rotmat = R_x(Const::pi / 4.0);

    Matrix expected(3, 3);
    expected(1,1) = 1.0; expected(1,2) = 0.0; expected(1,3) = 0.0;
    expected(2,1) = 0.0; expected(2,2) = 0.7071067811865475; expected(2,3) = 0.7071067811865475;
    expected(3,1) = 0.0; expected(3,2) = -0.7071067811865475; expected(3,3) = 0.7071067811865475;

    _assert(m_equals(rotmat, expected, 1e-6));
    return 0;
}

// Test para R_y
int m_r_y_01() {
    Matrix result = R_y(Const::pi / 4.0);

    Matrix expected(3, 3);
    expected(1,1) = 0.707106781186548; expected(1,2) = 0.0; expected(1,3) = -0.707106781186548;
    expected(2,1) = 0.0;               expected(2,2) = 1.0; expected(2,3) = 0.0;
    expected(3,1) = 0.707106781186548; expected(3,2) = 0.0; expected(3,3) = 0.707106781186548;

    _assert(m_equals(result, expected, 1e-12));
    return 0;
}

// Test para R_z
int m_r_z_01() {
    Matrix result = R_z(Const::pi / 4.0);

    Matrix expected(3, 3);
    expected(1, 1) = 0.707106781186548; expected(1, 2) = 0.707106781186548; expected(1, 3) = 0.0;
    expected(2, 1) = -0.707106781186548; expected(2, 2) = 0.707106781186548; expected(2, 3) = 0.0;
    expected(3, 1) = 0.0; expected(3, 2) = 0.0; expected(3, 3) = 1.0;


    _assert(m_equals(result, expected, 1e-9));

    return 0;
}

int m_sign_01() {
    double result = sign_(-5.0, 3.0);
    double expected = 5.0;
    _assert(fabs(result - expected) < 1e-10);
    return 0;
}

int m_timediff_01() {
    double UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC;
    timediff(-0.5, 37.0, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);
    double expected = 69.184;
    _assert(fabs(TT_UTC - expected) < 1e-6);
    return 0;
}

int m_azelpa_01() {
    Matrix s(3, 1);
    s(1,1) = 1.0; s(2,1) = 2.0; s(3,1) = 3.0;
    double Az, El;
    Matrix dAds(1, 3), dEds(1, 3);
    AzElPa(s, Az, El, dAds, dEds);
    double expected = 0.4636476090008061;
    _assert(fabs(Az - expected) < 1e-6);
    return 0;
}

int m_iers_01() {
	
    double Mjd_UTC = 49746.1163541665;      //tengo que cambiar por variable
    char interp = 'l';


    double x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC;

    // Llamar a IERS
    IERS(Mjd_UTC, interp, x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC);


    double expected_x_pole = -5.59378724204105e-07;  // rad
    double expected_y_pole = 2.33559834147171e-06;  // rad
    double expected_UT1_UTC = 0.325747632958789;    // s
    double expected_LOD = 0.0027269897187419;       // s
    double expected_dpsi = -1.16882953161739e-07;   // rad
    double expected_deps = -2.47835061986699e-08;   // rad
    double expected_dx_pole = -8.43027359620098e-10; // rad
    double expected_dy_pole = -1.56811369104134e-09; // rad
    double expected_TAI_UTC = 29.0;                 // s


    _assert(fabs(expected_x_pole - x_pole) < 1e-15);
    _assert(fabs(expected_y_pole - y_pole) < 1e-15);
    _assert(fabs(expected_UT1_UTC - UT1_UTC) < 1e-10);
    _assert(fabs(expected_LOD - LOD) < 1e-10);
    _assert(fabs(expected_dpsi - dpsi) < 1e-15);
    _assert(fabs(expected_deps - deps) < 1e-15);
    _assert(fabs(expected_dx_pole - dx_pole) < 1e-20);
    _assert(fabs(expected_dy_pole - dy_pole) < 1e-15);
    _assert(fabs(expected_TAI_UTC - TAI_UTC) < 1e-10);

    return 0;
}

int m_legendre_01() {
    int n = 3;
    int m = 3;
    double fi = Const::pi/4.0; // pi/4
    Matrix pnm(n+1, m+1);
    Matrix dpnm(n+1, m+1);
    Legendre(n, m, fi, pnm, dpnm);
    double expected = 1.224744871391589;
    _assert(fabs(pnm(2,2) - expected) < 1e-6);
    return 0;
}

int m_nutangles_01() {
    double Mjd_TT = 60355.0;
    double dpsi, deps;
    NutAngles(Mjd_TT, dpsi, deps);
    double expected = -2.022464e-05;
    _assert(fabs(dpsi - expected) < 1e-6);
    return 0;
}


int m_timeupdate_01() {
    Matrix P(2, 2);
    P(1,1) = 1.0; P(1,2) = 0.5; P(2,1) = 0.5; P(2,2) = 2.0;
    Matrix Phi(2, 2);
    Phi(1,1) = 1.0; Phi(1,2) = 0.0; Phi(2,1) = 0.0; Phi(2,2) = 1.0;
    Matrix Qdt(2, 2);
    Qdt(1,1) = 0.1; Qdt(1,2) = 0.0; Qdt(2,1) = 0.0; Qdt(2,2) = 0.2;
    Matrix P_new = TimeUpdate(P, Phi, Qdt);
    double expected = 1.1;
    _assert(fabs(P_new(1,1) - expected) < 1e-6);
    return 0;
}

int m_accel_harmonic_01() {
	eop19620101(21413);
	GGM03S(181);
	DE430Coeff(2285, 1020);

    Matrix r(3, 1);
    r(1, 1) = 7000e3;  // 7000 km = 7,000,000 m
    r(2, 1) = 0.0;
    r(3, 1) = 0.0;


    Matrix E = eye(3);


    int n_max = 2;
    int m_max = 2;


    Matrix R = AccelHarmonic(r, E, n_max, m_max);


    Matrix expected(3, 1);
    expected(1, 1) = -8.14576607065686;      // m/s²
    expected(2, 1) = -3.66267894892037e-05;  // m/s²
    expected(3, 1) = -5.84508413583961e-09;  // m/s²


    _assert(m_equals(expected, R, 1e-10));
    return 0;
}


// Test para EqnEquinox
int m_eqn_equinox_01() {
    // Inicializar datos globales si es necesario
    eop19620101(21413); // Asegurar que eopdata esté inicializado

    // Definir Mjd_TT
    double Mjd_TT = 49746.1163541665;

    // Calcular EqE
    double EqE = EqnEquinox(Mjd_TT);

    // Valor esperado de MATLAB
    double expected_EqE = 5.716767215940e-05; // rad

    // Verificar con tolerancia
    _assert(fabs(EqE - expected_EqE) < 1e-10);

    return 0;
}



// Test for JPL_Eph_DE430
int m_jpl_eph_de430_01() {
    // Initialize global PC matrix
    DE430Coeff(1, 1020); // Initialize PC with 1 row, 1020 columns
    PC(1, 1) = 2451545.0; // JD_start
    PC(1, 2) = 2451577.0; // JD_end

    // Earth coefficients (231:308)
    for (int i = 0; i < 13; ++i) {
        PC(1, 231 + i) = i + 1; // Cx_Earth, subinterval 1
        PC(1, 244 + i) = 14 + i; // Cy_Earth
        PC(1, 257 + i) = 27 + i; // Cz_Earth
        PC(1, 270 + i) = 40 + i; // Cx_Earth, subinterval 2
        PC(1, 283 + i) = 53 + i; // Cy_Earth
        PC(1, 296 + i) = 66 + i; // Cz_Earth
    }

    // Sun coefficients (753:818)
    for (int i = 0; i < 11; ++i) {
        PC(1, 753 + i) = i + 1; // Cx_Sun, subinterval 1
        PC(1, 764 + i) = 12 + i; // Cy_Sun
        PC(1, 775 + i) = 23 + i; // Cz_Sun
        PC(1, 786 + i) = 34 + i; // Cx_Sun, subinterval 2
        PC(1, 797 + i) = 45 + i; // Cy_Sun
        PC(1, 808 + i) = 56 + i; // Cz_Sun
    }

    // Moon coefficients (441:752, 8 subintervals)
    for (int s = 0; s < 8; ++s) {
        int base = 441 + s * 39;
        for (int i = 0; i < 13; ++i) {
            PC(1, base + i) = i + 1; // Cx_Moon
            PC(1, base + 13 + i) = 14 + i; // Cy_Moon
            PC(1, base + 26 + i) = 27 + i; // Cz_Moon
        }
    }

    // Mercury coefficients (3:45, subinterval 1 only)
    for (int i = 0; i < 14; ++i) {
        PC(1, 3 + i) = i + 1; // Cx_Mercury
        PC(1, 17 + i) = 15 + i; // Cy_Mercury
        PC(1, 31 + i) = 29 + i; // Cz_Mercury
    }

    // Test Mjd_TDB values
    double Mjd_TDB_values[] = {51544.5, 51560.5, 51548.5};
    for (int idx = 0; idx < 3; ++idx) {
        double Mjd_TDB = Mjd_TDB_values[idx];

        // Initialize output matrices
        Matrix r_Mercury(3, 1), r_Venus(3, 1), r_Earth(3, 1), r_Mars(3, 1),
               r_Jupiter(3, 1), r_Saturn(3, 1), r_Uranus(3, 1), r_Neptune(3, 1),
               r_Pluto(3, 1), r_Moon(3, 1), r_Sun(3, 1);

        // Call JPL_Eph_DE430
        JPL_Eph_DE430(Mjd_TDB, r_Mercury, r_Venus, r_Earth, r_Mars, r_Jupiter,
                      r_Saturn, r_Uranus, r_Neptune, r_Pluto, r_Moon, r_Sun);

        // Expected values
        Matrix exp_r_Moon(3, 1), exp_r_Earth(3, 1), exp_r_Sun(3, 1), exp_r_Mercury(3, 1),
               exp_r_Venus(3, 1), exp_r_Mars(3, 1), exp_r_Jupiter(3, 1), exp_r_Saturn(3, 1),
               exp_r_Uranus(3, 1), exp_r_Neptune(3, 1), exp_r_Pluto(3, 1);
        double EMRAT1 = 1.0 / (1.0 + 81.30056907419062);
        if (idx == 0 || idx == 2) { // dt=0 or dt=4
            exp_r_Moon(1, 1) = 1000; exp_r_Moon(2, 1) = 14000; exp_r_Moon(3, 1) = 27000;
            exp_r_Earth = exp_r_Moon - exp_r_Moon * EMRAT1;
            exp_r_Sun(1, 1) = 1000 - exp_r_Earth(1, 1);
            exp_r_Sun(2, 1) = 12000 - exp_r_Earth(2, 1);
            exp_r_Sun(3, 1) = 23000 - exp_r_Earth(3, 1);
            exp_r_Mercury(1, 1) = 1000 - exp_r_Earth(1, 1);
            exp_r_Mercury(2, 1) = 15000 - exp_r_Earth(2, 1);
            exp_r_Mercury(3, 1) = 29000 - exp_r_Earth(3, 1);
            // Others are zero
        } else { // dt=16
            exp_r_Moon(1, 1) = 1000; exp_r_Moon(2, 1) = 14000; exp_r_Moon(3, 1) = 27000;
            exp_r_Earth(1, 1) = 40000; exp_r_Earth(2, 1) = 53000; exp_r_Earth(3, 1) = 66000;
            exp_r_Earth = exp_r_Earth - exp_r_Moon * EMRAT1;
            exp_r_Sun(1, 1) = 34000 - exp_r_Earth(1, 1);
            exp_r_Sun(2, 1) = 45000 - exp_r_Earth(2, 1);
            exp_r_Sun(3, 1) = 56000 - exp_r_Earth(3, 1);
            exp_r_Mercury(1, 1) = 0 - exp_r_Earth(1, 1);
            exp_r_Mercury(2, 1) = 0 - exp_r_Earth(2, 1);
            exp_r_Mercury(3, 1) = 0 - exp_r_Earth(3, 1);
            // Others are zero
        }
		
			        std::cout << "AQUI LLEGAStest";

        // Verify results with print statements
        printf("Mjd_TDB = %.1f\n", Mjd_TDB);
        printf("r_Moon: [%.6f, %.6f, %.6f]\n", r_Moon(1, 1), r_Moon(2, 1), r_Moon(3, 1));
        printf("exp_r_Moon: [%.6f, %.6f, %.6f]\n", exp_r_Moon(1, 1), exp_r_Moon(2, 1), exp_r_Moon(3, 1));
        _assert(m_equals(r_Moon, exp_r_Moon, 1e-6));

        printf("r_Earth: [%.6f, %.6f, %.6f]\n", r_Earth(1, 1), r_Earth(2, 1), r_Earth(3, 1));
        printf("exp_r_Earth: [%.6f, %.6f, %.6f]\n", exp_r_Earth(1, 1), exp_r_Earth(2, 1), exp_r_Earth(3, 1));
        _assert(m_equals(r_Earth, exp_r_Earth, 1e-6));

        printf("r_Sun: [%.6f, %.6f, %.6f]\n", r_Sun(1, 1), r_Sun(2, 1), r_Sun(3, 1));
        printf("exp_r_Sun: [%.6f, %.6f, %.6f]\n", exp_r_Sun(1, 1), exp_r_Sun(2, 1), exp_r_Sun(3, 1));
        _assert(m_equals(r_Sun, exp_r_Sun, 1e-6));

        printf("r_Mercury: [%.6f, %.6f, %.6f]\n", r_Mercury(1, 1), r_Mercury(2, 1), r_Mercury(3, 1));
        printf("exp_r_Mercury: [%.6f, %.6f, %.6f]\n", exp_r_Mercury(1, 1), exp_r_Mercury(2, 1), exp_r_Mercury(3, 1));
        _assert(m_equals(r_Mercury, exp_r_Mercury, 1e-6));

        printf("r_Venus: [%.6f, %.6f, %.6f]\n", r_Venus(1, 1), r_Venus(2, 1), r_Venus(3, 1));
        printf("exp_r_Venus: [%.6f, %.6f, %.6f]\n", exp_r_Venus(1, 1), exp_r_Venus(2, 1), exp_r_Venus(3, 1));
        _assert(m_equals(r_Venus, exp_r_Venus, 1e-6));

        printf("r_Mars: [%.6f, %.6f, %.6f]\n", r_Mars(1, 1), r_Mars(2, 1), r_Mars(3, 1));
        printf("exp_r_Mars: [%.6f, %.6f, %.6f]\n", exp_r_Mars(1, 1), exp_r_Mars(2, 1), exp_r_Mars(3, 1));
        _assert(m_equals(r_Mars, exp_r_Mars, 1e-6));

        printf("r_Jupiter: [%.6f, %.6f, %.6f]\n", r_Jupiter(1, 1), r_Jupiter(2, 1), r_Jupiter(3, 1));
        printf("exp_r_Jupiter: [%.6f, %.6f, %.6f]\n", exp_r_Jupiter(1, 1), exp_r_Jupiter(2, 1), exp_r_Jupiter(3, 1));
        _assert(m_equals(r_Jupiter, exp_r_Jupiter, 1e-6));

        printf("r_Saturn: [%.6f, %.6f, %.6f]\n", r_Saturn(1, 1), r_Saturn(2, 1), r_Saturn(3, 1));
        printf("exp_r_Saturn: [%.6f, %.6f, %.6f]\n", exp_r_Saturn(1, 1), exp_r_Saturn(2, 1), exp_r_Saturn(3, 1));
        _assert(m_equals(r_Saturn, exp_r_Saturn, 1e-6));

        printf("r_Uranus: [%.6f, %.6f, %.6f]\n", r_Uranus(1, 1), r_Uranus(2, 1), r_Uranus(3, 1));
        printf("exp_r_Uranus: [%.6f, %.6f, %.6f]\n", exp_r_Uranus(1, 1), exp_r_Uranus(2, 1), exp_r_Uranus(3, 1));
        _assert(m_equals(r_Uranus, exp_r_Uranus, 1e-6));

        printf("r_Neptune: [%.6f, %.6f, %.6f]\n", r_Neptune(1, 1), r_Neptune(2, 1), r_Neptune(3, 1));
        printf("exp_r_Neptune: [%.6f, %.6f, %.6f]\n", exp_r_Neptune(1, 1), exp_r_Neptune(2, 1), exp_r_Neptune(3, 1));
        _assert(m_equals(r_Neptune, exp_r_Neptune, 1e-6));

        printf("r_Pluto: [%.6f, %.6f, %.6f]\n", r_Pluto(1, 1), r_Pluto(2, 1), r_Pluto(3, 1));
        printf("exp_r_Pluto: [%.6f, %.6f, %.6f]\n", exp_r_Pluto(1, 1), exp_r_Pluto(2, 1), exp_r_Pluto(3, 1));
        _assert(m_equals(r_Pluto, exp_r_Pluto, 1e-6));
    }

    return 0;
}














int all_tests() {
    _verify(m_constructor_01);
    _verify(m_access_01);
    _verify(m_sum_01);
    _verify(m_sub_01);
    _verify(m_mult_01);
    _verify(m_div_01);
    _verify(m_assign_01);
    // No test para operator << (se prueba implícitamente)
    _verify(m_zeros_01);
    _verify(m_eye_01);
    _verify(m_transpose_01);
    _verify(m_inv_01);
    _verify(m_scalar_add_01);
    _verify(m_scalar_sub_01);
    _verify(m_scalar_mult_01);
    _verify(m_scalar_div_01);
    _verify(m_constructor_n_01);
    _verify(m_vector_access_01);
    _verify(m_zeros_n_01);
    _verify(m_norm_01);
    _verify(m_dot_01);
    _verify(m_cross_01);
    _verify(m_extract_vector_01);
    _verify(m_union_vector_01);
    _verify(m_extract_row_01);
    _verify(m_extract_column_01);
    _verify(m_assign_row_01);
    _verify(m_assign_column_01);
	_verify(m_accel_point_mass_01);
	_verify(m_cheb3d_01);             
	_verify(m_ecc_anom_01);
	_verify(m_frac_01);
	_verify(m_mean_obliquity_01);       //test num 31
	_verify(m_mjday_01);
	_verify(m_mjday_tdb_01);   //33
	_verify(m_position_01);      //34    
	_verify(m_r_x_01);
	_verify(m_r_y_01);
	_verify(m_r_z_01);           //37
	// No test para SAT_Const (solo son variables)
	_verify(m_sign_01);
	_verify(m_timediff_01);
	_verify(m_azelpa_01);        //40
	_verify(m_iers_01);
	_verify(m_legendre_01);          //42
	_verify(m_nutangles_01);
	_verify(m_timeupdate_01);         //44+1
	_verify(m_accel_harmonic_01);     //46
	_verify(m_eqn_equinox_01);        //47
	_verify(m_jpl_eph_de430_01);
	
	
	
	

    return 0;
}

int main() {
	
	
    eop19620101(21413);
    GGM03S(181);
    DE430Coeff(2285, 1020);
	
	
	
    int result = all_tests();

    if (result == 0)
        printf("PASSED\n");

    printf("Tests run: %d\n", tests_run);

    return result != 0;
}