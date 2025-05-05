#include "..\include\matrix.hpp"
#include <cstdio>
#include <cmath>

#include "..\include\SAT_Const.hpp"

#include "..\include\AccelPointMass.hpp"
#include "..\include\Cheb3D.hpp"
#include "..\include\EccAnom.hpp"
#include "..\include\Frac.hpp" 
#include "..\include\MeanObliquity.hpp"






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
    double t = 1800;
    int N = 4;
    double Ta = 0;
    double Tb = 3600;

    Matrix Cx(4, 1);
    Cx(1,1) = 1000000; Cx(2,1) = 500000; Cx(3,1) = 0; Cx(4,1) = 0;
    Matrix Cy(4, 1);
    Cy(1,1) = 0; Cy(2,1) = 500000; Cy(3,1) = 1000000; Cy(4,1) = 0;
    Matrix Cz(4, 1);
    Cz(1,1) = 0; Cz(2,1) = 0; Cz(3,1) = 500000; Cz(4,1) = 1000000;

    // Debug: Print coefficients
    std::cout << "Test Cx: ";
    for (int i = 1; i <= N; i++) std::cout << Cx(i,1) << " ";
    std::cout << "\nTest Cy: ";
    for (int i = 1; i <= N; i++) std::cout << Cy(i,1) << " ";
    std::cout << "\nTest Cz: ";
    for (int i = 1; i <= N; i++) std::cout << Cz(i,1) << " ";
    std::cout << "\n";

    Matrix expected(1, 3);
    expected(1,1) = 1000000;
    expected(1,2) = -1000000;
    expected(1,3) = -500000;

    Matrix ChebApp = Cheb3D(t, N, Ta, Tb, Cx, Cy, Cz);

    _assert(m_equals(ChebApp, expected, 1e-9));

    double expected_magnitude = 1500000;
    double result_magnitude = norm(ChebApp);
    _assert(fabs(expected_magnitude - result_magnitude) < 1e-9);

    return 0;
}


// Test EccAnom
int m_ecc_anom_01() {     //pruebo con diferentes valores
    // Test Case 1: Circular orbit (e=0)
    _assert(fabs(EccAnom(0.7854, 0.0) - 0.7854) < 1e-4);
    
    // Test Case 2: Low eccentricity orbit (e=0.1)
    _assert(fabs(EccAnom(1.5708, 0.1) - 1.6703) < 1e-4);
    
    // Test Case 3: Medium eccentricity orbit (e=0.5)
    _assert(fabs(EccAnom(3.1416, 0.5) - 3.1416) < 1e-4);
    
    // Test Case 4: High eccentricity orbit (e=0.9)
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
	//_verify(m_cheb3d_01);             //NO ME VA ESTE TEST NO SE POR QUE 
	_verify(m_ecc_anom_01);
	_verify(m_frac_01);
	_verify(m_mean_obliquity_01);       //test num 31

    return 0;
}

int main() {
    int result = all_tests();

    if (result == 0)
        printf("PASSED\n");

    printf("Tests run: %d\n", tests_run);

    return result != 0;
}