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
    double Mjd_TDB = 60348.0;

    // Initialize output matrices
    Matrix r_Mercury(1, 3), r_Venus(1, 3), r_Earth(1, 3), r_Mars(1, 3),
           r_Jupiter(1, 3), r_Saturn(1, 3), r_Uranus(1, 3), r_Neptune(1, 3),
           r_Pluto(1, 3), r_Moon(1, 3), r_Sun(1, 3);

    // Call JPL_Eph_DE430
    JPL_Eph_DE430(Mjd_TDB, r_Mercury, r_Venus, r_Earth, r_Mars, r_Jupiter,
                  r_Saturn, r_Uranus, r_Neptune, r_Pluto, r_Moon, r_Sun);

    // Expected values
    Matrix expected_r_Mercury(1, 3);
    expected_r_Mercury(1, 1) = 112623958311.779;
    expected_r_Mercury(1, 2) = -150868623524.381;
    expected_r_Mercury(1, 3) = -71779751443.9542;
    _assert(m_equals(r_Mercury, expected_r_Mercury, 1e-3));
    std::cout << "Test passed: r_Mercury\n";

    Matrix expected_r_Venus(1, 3);
    expected_r_Venus(1, 1) = 67979665801.0159;
    expected_r_Venus(1, 2) = -182040469650.972;
    expected_r_Venus(1, 3) = -77754531631.0774;
    _assert(m_equals(r_Venus, expected_r_Venus, 1e-3));
    std::cout << "Test passed: r_Venus\n";

    Matrix expected_r_Earth(1, 3);
    expected_r_Earth(1, 1) = -111448106096.032;
    expected_r_Earth(1, 2) = 89490608943.6678;
    expected_r_Earth(1, 3) = 38828100871.8833;
    _assert(m_equals(r_Earth, expected_r_Earth, 1e-3));
    std::cout << "Test passed: r_Earth\n";

    Matrix expected_r_Mars(1, 3);
    expected_r_Mars(1, 1) = 148466860011.68;
    expected_r_Mars(1, 2) = -281603620116.961;
    expected_r_Mars(1, 3) = -127931017643.326;
    _assert(m_equals(r_Mars, expected_r_Mars, 1e-3));
    std::cout << "Test passed: r_Mars\n";

    Matrix expected_r_Jupiter(1, 3);
    expected_r_Jupiter(1, 1) = 600849652579.016;
    expected_r_Jupiter(1, 2) = 431907465367.823;
    expected_r_Jupiter(1, 3) = 172747882009.303;
    _assert(m_equals(r_Jupiter, expected_r_Jupiter, 1e-2));
    std::cout << "Test passed: r_Jupiter\n";

    Matrix expected_r_Saturn(1, 3);
    expected_r_Saturn(1, 1) = 1466091955279.12;
    expected_r_Saturn(1, 2) = -555189637093.051;
    expected_r_Saturn(1, 3) = -289531886725.713;
    _assert(m_equals(r_Saturn, expected_r_Saturn, 1e-3));
    std::cout << "Test passed: r_Saturn\n";

    Matrix expected_r_Uranus(1, 3);
    expected_r_Uranus(1, 1) = 1928309175333.92;
    expected_r_Uranus(1, 2) = 2027905259796.94;
    expected_r_Uranus(1, 3) = 862836636210.81;
    _assert(m_equals(r_Uranus, expected_r_Uranus, 1e-2));
    std::cout << "Test passed: r_Uranus\n";

    Matrix expected_r_Neptune(1, 3);
    expected_r_Neptune(1, 1) = 4575619355154.4;
    expected_r_Neptune(1, 2) = -280382588638.193;
    expected_r_Neptune(1, 3) = -228103255917.356;
    _assert(m_equals(r_Neptune, expected_r_Neptune, 1e-2));
    std::cout << "Test passed: r_Neptune\n";

    Matrix expected_r_Pluto(1, 3);
    expected_r_Pluto(1, 1) = 2700803532076.77;
    expected_r_Pluto(1, 2) = -4144525986799.46;
    expected_r_Pluto(1, 3) = -2084446068792.87;
    _assert(m_equals(r_Pluto, expected_r_Pluto, 1e-1));
    std::cout << "Test passed: r_Pluto\n";

    Matrix expected_r_Moon(1, 3);
    expected_r_Moon(1, 1) = 130639413.73261;
    expected_r_Moon(1, 2) = -298652884.800457;
    expected_r_Moon(1, 3) = -164607636.963072;
    _assert(m_equals(r_Moon, expected_r_Moon, 1e-4));
    std::cout << "Test passed: r_Moon\n";

    Matrix expected_r_Sun(1, 3);
    expected_r_Sun(1, 1) = 110284675128.512;
    expected_r_Sun(1, 2) = -89937859346.5398;
    expected_r_Sun(1, 3) = -38988031316.0913;
    _assert(m_equals(r_Sun, expected_r_Sun, 1e-3));
    std::cout << "Test passed: r_Sun\n";

    std::cout << "All tests passed successfully!\n";
    return 0;
}

// Test para EqnEquinox
int m_LTC_01() {
    // Test inputs
    double lon = Const::pi / 4.0; // 45 degrees
    double lat = Const::pi / 6.0; // 30 degrees

    // Call LTC
    Matrix M = LTC(lon, lat);

    // Expected transformation matrix
    Matrix expected_M(3, 3);
    expected_M(1, 1) = -0.7071067811865475;  expected_M(1, 2) = 0.7071067811865475;  expected_M(1, 3) = 0.0;
    expected_M(2, 1) = -0.3535533905932737;  expected_M(2, 2) = -0.3535533905932737; expected_M(2, 3) = 0.8660254037844386;
    expected_M(3, 1) = 0.6123724356957945;  expected_M(3, 2) = 0.6123724356957945; expected_M(3, 3) = 0.5;

    // Verify result
    _assert(m_equals(M, expected_M, 1e-6));

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
	_verify(m_jpl_eph_de430_01);          //48
	_verify(m_LTC_01);                           //  49 por si da error revisar este test
	
	
	
	
	

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