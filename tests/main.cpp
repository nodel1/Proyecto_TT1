#include "..\include\global.hpp"
#include "..\include\matrix.hpp"
#include <iostream>

using namespace std;

int main() {
	
	eop19620101(4);   //c =214123
	cout <<"eopdata\n" << eopdata << "\n";
	
	
	
	
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
	
	
	
    Matrix M1(3, 2);
	M1(1,1) = 5;
	
    Matrix M2(3, 2);
	M2(1,1) = -3;
	
    Matrix M3 = M1 - M2;

    cout << "M1\n" << M1 << "\n";
    cout << "M2\n" << M2 << "\n";
    cout << "M3\n" << M3 << "\n";
	
	cout << M1(1,1) << "\n";

    return 0;
}