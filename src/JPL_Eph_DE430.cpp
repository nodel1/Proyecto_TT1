// JPL_Eph_DE430.cpp
#include "..\include\JPL_Eph_DE430.hpp"
#include "..\include\Cheb3D.hpp"
#include "..\include\matrix.hpp"
#include "..\include\global.hpp"
#include "..\include\SAT_Const.hpp"

void JPL_Eph_DE430(double Mjd_TDB, Matrix& r_Mercury, Matrix& r_Venus, Matrix& r_Earth,
                   Matrix& r_Mars, Matrix& r_Jupiter, Matrix& r_Saturn, Matrix& r_Uranus,
                   Matrix& r_Neptune, Matrix& r_Pluto, Matrix& r_Moon, Matrix& r_Sun) {
    // Julian Date
    double JD = Mjd_TDB + 2400000.5;

        std::cout << "AQUI LLEGAS1";
		
    // Find interval in PC
    int i = 0;
    for (int j = 1; j <= PC.n_row; j++) {
        if (PC(j, 1) <= JD && JD <= PC(j, 2)) {
            i = j;
            break;
        }
    }
    if (i == 0) {
        // Error: interval not found
        return;
    }


        std::cout << "AQUI LLEGAS2";
		
		
    // Extract row from PC
    Matrix PCtemp = extract_row(PC, i);

    // Time relative to start of interval
    double t1 = PCtemp(1, 1) - 2400000.5; // MJD at start
    double dt = Mjd_TDB - t1;

    // Earth
    int temp_aux = (270 - 231) / 13 + 1;
    Matrix temp = zeros(temp_aux, 1);
    for (int k = 1; k <= temp_aux; k++) {
        temp(k, 1) = (k - 1) * 13 + 231;
    }
    Matrix Cx_Earth(temp(2, 1) - temp(1, 1), 1);
    Matrix Cy_Earth(temp(3, 1) - temp(2, 1), 1);
    Matrix Cz_Earth(temp(4, 1) - temp(3, 1), 1);
    for (int col = temp(1, 1), idx = 1; col < temp(2, 1); col++, idx++) {
        Cx_Earth(idx, 1) = PCtemp(1, col);
    }
    for (int col = temp(2, 1), idx = 1; col < temp(3, 1); col++, idx++) {
        Cy_Earth(idx, 1) = PCtemp(1, col);
    }
    for (int col = temp(3, 1), idx = 1; col < temp(4, 1); col++, idx++) {
        Cz_Earth(idx, 1) = PCtemp(1, col);
    }
	
	
	        std::cout << "AQUI LLEGAS3";
	
    temp = temp + 39;
    Matrix Cx(temp(2, 1) - temp(1, 1), 1);
    Matrix Cy(temp(3, 1) - temp(2, 1), 1);
    Matrix Cz(temp(4, 1) - temp(3, 1), 1);
    for (int col = temp(1, 1), idx = 1; col < temp(2, 1); col++, idx++) {
        Cx(idx, 1) = PCtemp(1, col);
    }
    for (int col = temp(2, 1), idx = 1; col < temp(3, 1); col++, idx++) {
        Cy(idx, 1) = PCtemp(1, col);
    }
    for (int col = temp(3, 1), idx = 1; col < temp(4, 1); col++, idx++) {
        Cz(idx, 1) = PCtemp(1, col);
    }
    Matrix Cx_Earth_full(Cx_Earth.n_row + Cx.n_row, 1);
    Matrix Cy_Earth_full(Cy_Earth.n_row + Cy.n_row, 1);
    Matrix Cz_Earth_full(Cz_Earth.n_row + Cz.n_row, 1);
    for (int k = 1; k <= Cx_Earth.n_row; k++) {
        Cx_Earth_full(k, 1) = Cx_Earth(k, 1);
        Cy_Earth_full(k, 1) = Cy_Earth(k, 1);
        Cz_Earth_full(k, 1) = Cz_Earth(k, 1);
    }
    for (int k = 1; k <= Cx.n_row; k++) {
        Cx_Earth_full(Cx_Earth.n_row + k, 1) = Cx(k, 1);
        Cy_Earth_full(Cy_Earth.n_row + k, 1) = Cy(k, 1);
        Cz_Earth_full(Cz_Earth.n_row + k, 1) = Cz(k, 1);
    }
	
	
	        std::cout << "AQUI LLEGAS4";
	
	
	
    int j;
    double Mjd0;
    if (0 <= dt && dt <= 16) {
        j = 0;
        Mjd0 = t1;
    } else {
        j = 1;
        Mjd0 = t1 + 16 * j;
    }
    Matrix Cx_Earth_sub(1, 13);
    Matrix Cy_Earth_sub(1, 13);
    Matrix Cz_Earth_sub(1, 13);
    for (int k = 1; k <= 13; k++) {
        Cx_Earth_sub(1, k) = Cx_Earth_full(13 * j + k, 1);
        Cy_Earth_sub(1, k) = Cy_Earth_full(13 * j + k, 1);
        Cz_Earth_sub(1, k) = Cz_Earth_full(13 * j + k, 1);
    }
    r_Earth = Cheb3D(Mjd_TDB, 13, Mjd0, Mjd0 + 16, Cx_Earth_sub, Cy_Earth_sub, Cz_Earth_sub) * 1e3;



        std::cout << "AQUI LLEGAS5";

    // Moon
    temp_aux = (480 - 441) / 13 + 1;
    temp = zeros(temp_aux, 1);
    for (int k = 1; k <= temp_aux; k++) {
        temp(k, 1) = (k - 1) * 13 + 441;
    }
    Matrix Cx_Moon(temp(2, 1) - temp(1, 1), 1);
    Matrix Cy_Moon(temp(3, 1) - temp(2, 1), 1);
    Matrix Cz_Moon(temp(4, 1) - temp(3, 1), 1);
    for (int col = temp(1, 1), idx = 1; col < temp(2, 1); col++, idx++) {
        Cx_Moon(idx, 1) = PCtemp(1, col);
    }
    for (int col = temp(2, 1), idx = 1; col < temp(3, 1); col++, idx++) {
        Cy_Moon(idx, 1) = PCtemp(1, col);
    }
    for (int col = temp(3, 1), idx = 1; col < temp(4, 1); col++, idx++) {
        Cz_Moon(idx, 1) = PCtemp(1, col);
    }
    Matrix Cx_Moon_full = Cx_Moon;
    Matrix Cy_Moon_full = Cy_Moon;
    Matrix Cz_Moon_full = Cz_Moon;
	
	        std::cout << "AQUI LLEGAS6";
	
	
    for (int k = 1; k <= 7; k++) {
        temp = temp + 39;
        Matrix Cx(temp(2, 1) - temp(1, 1), 1);
        Matrix Cy(temp(3, 1) - temp(2, 1), 1);
        Matrix Cz(temp(4, 1) - temp(3, 1), 1);
        for (int col = temp(1, 1), idx = 1; col < temp(2, 1); col++, idx++) {
            Cx(idx, 1) = PCtemp(1, col);
        }
        for (int col = temp(2, 1), idx = 1; col < temp(3, 1); col++, idx++) {
            Cy(idx, 1) = PCtemp(1, col);
        }
        for (int col = temp(3, 1), idx = 1; col < temp(4, 1); col++, idx++) {
            Cz(idx, 1) = PCtemp(1, col);
        }
        Matrix Cx_temp(Cx_Moon_full.n_row + Cx.n_row, 1);
        Matrix Cy_temp(Cy_Moon_full.n_row + Cy.n_row, 1);
        Matrix Cz_temp(Cz_Moon_full.n_row + Cz.n_row, 1);
        for (int m = 1; m <= Cx_Moon_full.n_row; m++) {
            Cx_temp(m, 1) = Cx_Moon_full(m, 1);
            Cy_temp(m, 1) = Cy_Moon_full(m, 1);
            Cz_temp(m, 1) = Cz_Moon_full(m, 1);
        }
        for (int m = 1; m <= Cx.n_row; m++) {
            Cx_temp(Cx_Moon_full.n_row + m, 1) = Cx(m, 1);
            Cy_temp(Cy_Moon_full.n_row + m, 1) = Cy(m, 1);
            Cz_temp(Cz_Moon_full.n_row + m, 1) = Cz(m, 1);
        }
        Cx_Moon_full = Cx_temp;
        Cy_Moon_full = Cy_temp;
        Cz_Moon_full = Cz_temp;
    }
	
	
	        std::cout << "AQUI LLEGAS7";
	
	
    if (0 <= dt && dt <= 4) j = 0;
    else if (4 < dt && dt <= 8) j = 1;
    else if (8 < dt && dt <= 12) j = 2;
    else if (12 < dt && dt <= 16) j = 3;
    else if (16 < dt && dt <= 20) j = 4;
    else if (20 < dt && dt <= 24) j = 5;
    else if (24 < dt && dt <= 28) j = 6;
    else j = 7;
    Mjd0 = t1 + 4 * j;
    Matrix Cx_Moon_sub(1, 13);
    Matrix Cy_Moon_sub(1, 13);
    Matrix Cz_Moon_sub(1, 13);
    for (int k = 1; k <= 13; k++) {
        Cx_Moon_sub(1, k) = Cx_Moon_full(13 * j + k, 1);
        Cy_Moon_sub(1, k) = Cy_Moon_full(13 * j + k, 1);
        Cz_Moon_sub(1, k) = Cz_Moon_full(13 * j + k, 1);
    }
    r_Moon = Cheb3D(Mjd_TDB, 13, Mjd0, Mjd0 + 4, Cx_Moon_sub, Cy_Moon_sub, Cz_Moon_sub) * 1e3;



        std::cout << "AQUI LLEGAS8";

    // Sun
    temp_aux = (786 - 753) / 11 + 1;
    temp = zeros(temp_aux, 1);
    for (int k = 1; k <= temp_aux; k++) {
        temp(k, 1) = (k - 1) * 11 + 753;
    }
    Matrix Cx_Sun(temp(2, 1) - temp(1, 1), 1);
    Matrix Cy_Sun(temp(3, 1) - temp(2, 1), 1);
    Matrix Cz_Sun(temp(4, 1) - temp(3, 1), 1);
    for (int col = temp(1, 1), idx = 1; col < temp(2, 1); col++, idx++) {
        Cx_Sun(idx, 1) = PCtemp(1, col);
    }
    for (int col = temp(2, 1), idx = 1; col < temp(3, 1); col++, idx++) {
        Cy_Sun(idx, 1) = PCtemp(1, col);
    }
    for (int col = temp(3, 1), idx = 1; col < temp(4, 1); col++, idx++) {
        Cz_Sun(idx, 1) = PCtemp(1, col);
    }
    temp = temp + 33;
    Matrix Cx_temp(temp(2, 1) - temp(1, 1), 1);
    Matrix Cy_temp(temp(3, 1) - temp(2, 1), 1);
    Matrix Cz_temp(temp(4, 1) - temp(3, 1), 1);
    for (int col = temp(1, 1), idx = 1; col < temp(2, 1); col++, idx++) {
        Cx_temp(idx, 1) = PCtemp(1, col);
    }
    for (int col = temp(2, 1), idx = 1; col < temp(3, 1); col++, idx++) {
        Cy_temp(idx, 1) = PCtemp(1, col);
    }
    for (int col = temp(3, 1), idx = 1; col < temp(4, 1); col++, idx++) {
        Cz_temp(idx, 1) = PCtemp(1, col);
    }
    Matrix Cx_Sun_full(Cx_Sun.n_row + Cx_temp.n_row, 1);
    Matrix Cy_Sun_full(Cy_Sun.n_row + Cy_temp.n_row, 1);
    Matrix Cz_Sun_full(Cz_Sun.n_row + Cz_temp.n_row, 1);
    for (int k = 1; k <= Cx_Sun.n_row; k++) {
        Cx_Sun_full(k, 1) = Cx_Sun(k, 1);
        Cy_Sun_full(k, 1) = Cy_Sun(k, 1);
        Cz_Sun_full(k, 1) = Cz_Sun(k, 1);
    }
    for (int k = 1; k <= Cx_temp.n_row; k++) {
        Cx_Sun_full(Cx_Sun.n_row + k, 1) = Cx_temp(k, 1);
        Cy_Sun_full(Cy_Sun.n_row + k, 1) = Cy_temp(k, 1);
        Cz_Sun_full(Cz_Sun.n_row + k, 1) = Cz_temp(k, 1);
    }
    if (0 <= dt && dt <= 16) {
        j = 0;
        Mjd0 = t1;
    } else {
        j = 1;
        Mjd0 = t1 + 16 * j;
    }
    Matrix Cx_Sun_sub(1, 11);
    Matrix Cy_Sun_sub(1, 11);
    Matrix Cz_Sun_sub(1, 11);
    for (int k = 1; k <= 11; k++) {
        Cx_Sun_sub(1, k) = Cx_Sun_full(11 * j + k, 1);
        Cy_Sun_sub(1, k) = Cy_Sun_full(11 * j + k, 1);
        Cz_Sun_sub(1, k) = Cz_Sun_full(11 * j + k, 1);
    }
    r_Sun = Cheb3D(Mjd_TDB, 11, Mjd0, Mjd0 + 16, Cx_Sun_sub, Cy_Sun_sub, Cz_Sun_sub) * 1e3;


	        std::cout << "AQUI LLEGAS sol";

    // Mercury
    temp_aux = (45 - 3) / 14 + 1;
    temp = zeros(temp_aux, 1);
    for (int k = 1; k <= temp_aux; k++) {
        temp(k, 1) = (k - 1) * 14 + 3;
    }
    Matrix Cx_Mercury(temp(2, 1) - temp(1, 1), 1);
    Matrix Cy_Mercury(temp(3, 1) - temp(2, 1), 1);
    Matrix Cz_Mercury(temp(4, 1) - temp(3, 1), 1);
    for (int col = temp(1, 1), idx = 1; col < temp(2, 1); col++, idx++) {
        Cx_Mercury(idx, 1) = PCtemp(1, col);
    }
    for (int col = temp(2, 1), idx = 1; col < temp(3, 1); col++, idx++) {
        Cy_Mercury(idx, 1) = PCtemp(1, col);
    }
    for (int col = temp(3, 1), idx = 1; col < temp(4, 1); col++, idx++) {
        Cz_Mercury(idx, 1) = PCtemp(1, col);
    }
    Matrix Cx_Mercury_full = Cx_Mercury;
    Matrix Cy_Mercury_full = Cy_Mercury;
    Matrix Cz_Mercury_full = Cz_Mercury;
    for (int k = 1; k <= 3; k++) {
        temp = temp + 42;
        Matrix Cx(temp(2, 1) - temp(1, 1), 1);
        Matrix Cy(temp(3, 1) - temp(2, 1), 1);
        Matrix Cz(temp(4, 1) - temp(3, 1), 1);
        for (int col = temp(1, 1), idx = 1; col < temp(2, 1); col++, idx++) {
            Cx(idx, 1) = PCtemp(1, col);
        }
        for (int col = temp(2, 1), idx = 1; col < temp(3, 1); col++, idx++) {
            Cy(idx, 1) = PCtemp(1, col);
        }
        for (int col = temp(3, 1), idx = 1; col < temp(4, 1); col++, idx++) {
            Cz(idx, 1) = PCtemp(1, col);
        }
        Matrix Cx_temp(Cx_Mercury_full.n_row + Cx.n_row, 1);
        Matrix Cy_temp(Cy_Mercury_full.n_row + Cy.n_row, 1);
        Matrix Cz_temp(Cz_Mercury_full.n_row + Cz.n_row, 1);
        for (int m = 1; m <= Cx_Mercury_full.n_row; m++) {
            Cx_temp(m, 1) = Cx_Mercury_full(m, 1);
            Cy_temp(m, 1) = Cy_Mercury_full(m, 1);
            Cz_temp(m, 1) = Cz_Mercury_full(m, 1);
        }
        for (int m = 1; m <= Cx.n_row; m++) {
            Cx_temp(Cx_Mercury_full.n_row + m, 1) = Cx(m, 1);
            Cy_temp(Cy_Mercury_full.n_row + m, 1) = Cy(m, 1);
            Cz_temp(Cz_Mercury_full.n_row + m, 1) = Cz(m, 1);
        }
        Cx_Mercury_full = Cx_temp;
        Cy_Mercury_full = Cy_temp;
        Cz_Mercury_full = Cz_temp;
    }
    if (0 <= dt && dt <= 8) j = 0;
    else if (8 < dt && dt <= 16) j = 1;
    else if (16 < dt && dt <= 24) j = 2;
    else j = 3;
    Mjd0 = t1 + 8 * j;
    Matrix Cx_Mercury_sub(1, 14);
    Matrix Cy_Mercury_sub(1, 14);
    Matrix Cz_Mercury_sub(1, 14);
    for (int k = 1; k <= 14; k++) {
        Cx_Mercury_sub(1, k) = Cx_Mercury_full(14 * j + k, 1);
        Cy_Mercury_sub(1, k) = Cy_Mercury_full(14 * j + k, 1);
        Cz_Mercury_sub(1, k) = Cz_Mercury_full(14 * j + k, 1);
    }
    r_Mercury = Cheb3D(Mjd_TDB, 14, Mjd0, Mjd0 + 8, Cx_Mercury_sub, Cy_Mercury_sub, Cz_Mercury_sub) * 1e3;



	        std::cout << "AQUI LLEGASmercurio";
			
			
    // Venus
    temp_aux = (201 - 171) / 10 + 1;
    temp = zeros(temp_aux, 1);
    for (int k = 1; k <= temp_aux; k++) {
        temp(k, 1) = (k - 1) * 10 + 171;
    }
    Matrix Cx_Venus(temp(2, 1) - temp(1, 1), 1);
    Matrix Cy_Venus(temp(3, 1) - temp(2, 1), 1);
    Matrix Cz_Venus(temp(4, 1) - temp(3, 1), 1);
    for (int col = temp(1, 1), idx = 1; col < temp(2, 1); col++, idx++) {
        Cx_Venus(idx, 1) = PCtemp(1, col);
    }
    for (int col = temp(2, 1), idx = 1; col < temp(3, 1); col++, idx++) {
        Cy_Venus(idx, 1) = PCtemp(1, col);
    }
    for (int col = temp(3, 1), idx = 1; col < temp(4, 1); col++, idx++) {
        Cz_Venus(idx, 1) = PCtemp(1, col);
    }
    temp = temp + 30;
    Matrix Cx_temp_venus(temp(2, 1) - temp(1, 1), 1);
    Matrix Cy_temp_venus(temp(3, 1) - temp(2, 1), 1);
    Matrix Cz_temp_venus(temp(4, 1) - temp(3, 1), 1);
    for (int col = temp(1, 1), idx = 1; col < temp(2, 1); col++, idx++) {
        Cx_temp_venus(idx, 1) = PCtemp(1, col);
    }
    for (int col = temp(2, 1), idx = 1; col < temp(3, 1); col++, idx++) {
        Cy_temp_venus(idx, 1) = PCtemp(1, col);
    }
    for (int col = temp(3, 1), idx = 1; col < temp(4, 1); col++, idx++) {
        Cz_temp_venus(idx, 1) = PCtemp(1, col);
    }
    Matrix Cx_Venus_full(Cx_Venus.n_row + Cx_temp_venus.n_row, 1);
    Matrix Cy_Venus_full(Cy_Venus.n_row + Cy_temp_venus.n_row, 1);
    Matrix Cz_Venus_full(Cz_Venus.n_row + Cz_temp_venus.n_row, 1);
    for (int k = 1; k <= Cx_Venus.n_row; k++) {
        Cx_Venus_full(k, 1) = Cx_Venus(k, 1);
        Cy_Venus_full(k, 1) = Cy_Venus(k, 1);
        Cz_Venus_full(k, 1) = Cz_Venus(k, 1);
    }
    for (int k = 1; k <= Cx_temp_venus.n_row; k++) {
        Cx_Venus_full(Cx_Venus.n_row + k, 1) = Cx_temp_venus(k, 1);
        Cy_Venus_full(Cy_Venus.n_row + k, 1) = Cy_temp_venus(k, 1);
        Cz_Venus_full(Cz_Venus.n_row + k, 1) = Cz_temp_venus(k, 1);
    }
    if (0 <= dt && dt <= 16) {
        j = 0;
        Mjd0 = t1;
    } else {
        j = 1;
        Mjd0 = t1 + 16 * j;
    }
    Matrix Cx_Venus_sub(1, 10);
    Matrix Cy_Venus_sub(1, 10);
    Matrix Cz_Venus_sub(1, 10);
    for (int k = 1; k <= 10; k++) {
        Cx_Venus_sub(1, k) = Cx_Venus_full(10 * j + k, 1);
        Cy_Venus_sub(1, k) = Cy_Venus_full(10 * j + k, 1);
        Cz_Venus_sub(1, k) = Cz_Venus_full(10 * j + k, 1);
    }
    r_Venus = Cheb3D(Mjd_TDB, 10, Mjd0, Mjd0 + 16, Cx_Venus_sub, Cy_Venus_sub, Cz_Venus_sub) * 1e3;




	        std::cout << "AQUI LLEGASvenus";
			
			
			
			
    // Mars
    temp_aux = (342 - 309) / 11 + 1;
    temp = zeros(temp_aux, 1);
    for (int k = 1; k <= temp_aux; k++) {
        temp(k, 1) = (k - 1) * 11 + 309;
    }
    Matrix Cx_Mars(temp(2, 1) - temp(1, 1), 1);
    Matrix Cy_Mars(temp(3, 1) - temp(2, 1), 1);
    Matrix Cz_Mars(temp(4, 1) - temp(3, 1), 1);
    for (int col = temp(1, 1), idx = 1; col < temp(2, 1); col++, idx++) {
        Cx_Mars(idx, 1) = PCtemp(1, col);
    }
    for (int col = temp(2, 1), idx = 1; col < temp(3, 1); col++, idx++) {
        Cy_Mars(idx, 1) = PCtemp(1, col);
    }
    for (int col = temp(3, 1), idx = 1; col < temp(4, 1); col++, idx++) {
        Cz_Mars(idx, 1) = PCtemp(1, col);
    }
    j = 0;
    Mjd0 = t1;
    Matrix Cx_Mars_sub(1, 11);
    Matrix Cy_Mars_sub(1, 11);
    Matrix Cz_Mars_sub(1, 11);
    for (int k = 1; k <= 11; k++) {
        Cx_Mars_sub(1, k) = Cx_Mars(k, 1);
        Cy_Mars_sub(1, k) = Cy_Mars(k, 1);
        Cz_Mars_sub(1, k) = Cz_Mars(k, 1);
    }
    r_Mars = Cheb3D(Mjd_TDB, 11, Mjd0, Mjd0 + 32, Cx_Mars_sub, Cy_Mars_sub, Cz_Mars_sub) * 1e3;



	        std::cout << "AQUI LLEGASmarte";
			
			
			
    // Jupiter
    temp_aux = (366 - 342) / 8 + 1;
    temp = zeros(temp_aux, 1);
    for (int k = 1; k <= temp_aux; k++) {
        temp(k, 1) = (k - 1) * 8 + 342;
    }
    Matrix Cx_Jupiter(temp(2, 1) - temp(1, 1), 1);
    Matrix Cy_Jupiter(temp(3, 1) - temp(2, 1), 1);
    Matrix Cz_Jupiter(temp(4, 1) - temp(3, 1), 1);
    for (int col = temp(1, 1), idx = 1; col < temp(2, 1); col++, idx++) {
        Cx_Jupiter(idx, 1) = PCtemp(1, col);
    }
    for (int col = temp(2, 1), idx = 1; col < temp(3, 1); col++, idx++) {
        Cy_Jupiter(idx, 1) = PCtemp(1, col);
    }
    for (int col = temp(3, 1), idx = 1; col < temp(4, 1); col++, idx++) {
        Cz_Jupiter(idx, 1) = PCtemp(1, col);
    }
    j = 0;
    Mjd0 = t1;
    Matrix Cx_Jupiter_sub(1, 8);
    Matrix Cy_Jupiter_sub(1, 8);
    Matrix Cz_Jupiter_sub(1, 8);
    for (int k = 1; k <= 8; k++) {
        Cx_Jupiter_sub(1, k) = Cx_Jupiter(k, 1);
        Cy_Jupiter_sub(1, k) = Cy_Jupiter(k, 1);
        Cz_Jupiter_sub(1, k) = Cz_Jupiter(k, 1);
    }
    r_Jupiter = Cheb3D(Mjd_TDB, 8, Mjd0, Mjd0 + 32, Cx_Jupiter_sub, Cy_Jupiter_sub, Cz_Jupiter_sub) * 1e3;




	        std::cout << "AQUI LLEGASjupiter";
			
			
			
    // Saturn
    temp_aux = (387 - 366) / 7 + 1;
    temp = zeros(temp_aux, 1);
    for (int k = 1; k <= temp_aux; k++) {
        temp(k, 1) = (k - 1) * 7 + 366;
    }
    Matrix Cx_Saturn(temp(2, 1) - temp(1, 1), 1);
    Matrix Cy_Saturn(temp(3, 1) - temp(2, 1), 1);
    Matrix Cz_Saturn(temp(4, 1) - temp(3, 1), 1);
    for (int col = temp(1, 1), idx = 1; col < temp(2, 1); col++, idx++) {
        Cx_Saturn(idx, 1) = PCtemp(1, col);
    }
    for (int col = temp(2, 1), idx = 1; col < temp(3, 1); col++, idx++) {
        Cy_Saturn(idx, 1) = PCtemp(1, col);
    }
    for (int col = temp(3, 1), idx = 1; col < temp(4, 1); col++, idx++) {
        Cz_Saturn(idx, 1) = PCtemp(1, col);
    }
    j = 0;
    Mjd0 = t1;
    Matrix Cx_Saturn_sub(1, 7);
    Matrix Cy_Saturn_sub(1, 7);
    Matrix Cz_Saturn_sub(1, 7);
    for (int k = 1; k <= 7; k++) {
        Cx_Saturn_sub(1, k) = Cx_Saturn(k, 1);
        Cy_Saturn_sub(1, k) = Cy_Saturn(k, 1);
        Cz_Saturn_sub(1, k) = Cz_Saturn(k, 1);
    }
    r_Saturn = Cheb3D(Mjd_TDB, 7, Mjd0, Mjd0 + 32, Cx_Saturn_sub, Cy_Saturn_sub, Cz_Saturn_sub) * 1e3;

    // Uranus
    temp_aux = (405 - 387) / 6 + 1;
    temp = zeros(temp_aux, 1);
    for (int k = 1; k <= temp_aux; k++) {
        temp(k, 1) = (k - 1) * 6 + 387;
    }
    Matrix Cx_Uranus(temp(2, 1) - temp(1, 1), 1);
    Matrix Cy_Uranus(temp(3, 1) - temp(2, 1), 1);
    Matrix Cz_Uranus(temp(4, 1) - temp(3, 1), 1);
    for (int col = temp(1, 1), idx = 1; col < temp(2, 1); col++, idx++) {
        Cx_Uranus(idx, 1) = PCtemp(1, col);
    }
    for (int col = temp(2, 1), idx = 1; col < temp(3, 1); col++, idx++) {
        Cy_Uranus(idx, 1) = PCtemp(1, col);
    }
    for (int col = temp(3, 1), idx = 1; col < temp(4, 1); col++, idx++) {
        Cz_Uranus(idx, 1) = PCtemp(1, col);
    }
    j = 0;
    Mjd0 = t1;
    Matrix Cx_Uranus_sub(1, 6);
    Matrix Cy_Uranus_sub(1, 6);
    Matrix Cz_Uranus_sub(1, 6);
    for (int k = 1; k <= 6; k++) {
        Cx_Uranus_sub(1, k) = Cx_Uranus(k, 1);
        Cy_Uranus_sub(1, k) = Cy_Uranus(k, 1);
        Cz_Uranus_sub(1, k) = Cz_Uranus(k, 1);
    }
    r_Uranus = Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0 + 32, Cx_Uranus_sub, Cy_Uranus_sub, Cz_Uranus_sub) * 1e3;

    // Neptune
    temp_aux = (423 - 405) / 6 + 1;
    temp = zeros(temp_aux, 1);
    for (int k = 1; k <= temp_aux; k++) {
        temp(k, 1) = (k - 1) * 6 + 405;
    }
    Matrix Cx_Neptune(temp(2, 1) - temp(1, 1), 1);
    Matrix Cy_Neptune(temp(3, 1) - temp(2, 1), 1);
    Matrix Cz_Neptune(temp(4, 1) - temp(3, 1), 1);
    for (int col = temp(1, 1), idx = 1; col < temp(2, 1); col++, idx++) {
        Cx_Neptune(idx, 1) = PCtemp(1, col);
    }
    for (int col = temp(2, 1), idx = 1; col < temp(3, 1); col++, idx++) {
        Cy_Neptune(idx, 1) = PCtemp(1, col);
    }
    for (int col = temp(3, 1), idx = 1; col < temp(4, 1); col++, idx++) {
        Cz_Neptune(idx, 1) = PCtemp(1, col);
    }
    j = 0;
    Mjd0 = t1;
    Matrix Cx_Neptune_sub(1, 6);
    Matrix Cy_Neptune_sub(1, 6);
    Matrix Cz_Neptune_sub(1, 6);
    for (int k = 1; k <= 6; k++) {
        Cx_Neptune_sub(1, k) = Cx_Neptune(k, 1);
        Cy_Neptune_sub(1, k) = Cy_Neptune(k, 1);
        Cz_Neptune_sub(1, k) = Cz_Neptune(k, 1);
    }
    r_Neptune = Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0 + 32, Cx_Neptune_sub, Cy_Neptune_sub, Cz_Neptune_sub) * 1e3;

    // Pluto
    temp_aux = (441 - 423) / 6 + 1;
    temp = zeros(temp_aux, 1);
    for (int k = 1; k <= temp_aux; k++) {
        temp(k, 1) = (k - 1) * 6 + 423;
    }
    Matrix Cx_Pluto(temp(2, 1) - temp(1, 1), 1);
    Matrix Cy_Pluto(temp(3, 1) - temp(2, 1), 1);
    Matrix Cz_Pluto(temp(4, 1) - temp(3, 1), 1);
    for (int col = temp(1, 1), idx = 1; col < temp(2, 1); col++, idx++) {
        Cx_Pluto(idx, 1) = PCtemp(1, col);
    }
    for (int col = temp(2, 1), idx = 1; col < temp(3, 1); col++, idx++) {
        Cy_Pluto(idx, 1) = PCtemp(1, col);
    }
    for (int col = temp(3, 1), idx = 1; col < temp(4, 1); col++, idx++) {
        Cz_Pluto(idx, 1) = PCtemp(1, col);
    }
    j = 0;
    Mjd0 = t1;
    Matrix Cx_Pluto_sub(1, 6);
    Matrix Cy_Pluto_sub(1, 6);
    Matrix Cz_Pluto_sub(1, 6);
    for (int k = 1; k <= 6; k++) {
        Cx_Pluto_sub(1, k) = Cx_Pluto(k, 1);
        Cy_Pluto_sub(1, k) = Cy_Pluto(k, 1);
        Cz_Pluto_sub(1, k) = Cz_Pluto(k, 1);
    }
    r_Pluto = Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0 + 32, Cx_Pluto_sub, Cy_Pluto_sub, Cz_Pluto_sub) * 1e3;



	        std::cout << "AQUI LLEGASpluto";
			
			
			
    // Nutations
    temp_aux = (839 - 819) / 10 + 1;
    temp = zeros(temp_aux, 1);
    for (int k = 1; k <= temp_aux; k++) {
        temp(k, 1) = (k - 1) * 10 + 819;
    }
    Matrix Cx_Nutations(temp(2, 1) - temp(1, 1), 1);
    Matrix Cy_Nutations(temp(3, 1) - temp(2, 1), 1);
    for (int col = temp(1, 1), idx = 1; col < temp(2, 1); col++, idx++) {
        Cx_Nutations(idx, 1) = PCtemp(1, col);
    }
    for (int col = temp(2, 1), idx = 1; col < temp(3, 1); col++, idx++) {
        Cy_Nutations(idx, 1) = PCtemp(1, col);
    }
    Matrix Cx_Nutations_full = Cx_Nutations;
    Matrix Cy_Nutations_full = Cy_Nutations;
    for (int k = 1; k <= 3; k++) {
        temp = temp + 20;
        Matrix Cx(temp(2, 1) - temp(1, 1), 1);
        Matrix Cy(temp(3, 1) - temp(2, 1), 1);
        for (int col = temp(1, 1), idx = 1; col < temp(2, 1); col++, idx++) {
            Cx(idx, 1) = PCtemp(1, col);
        }
        for (int col = temp(2, 1), idx = 1; col < temp(3, 1); col++, idx++) {
            Cy(idx, 1) = PCtemp(1, col);
        }
        Matrix Cx_temp(Cx_Nutations_full.n_row + Cx.n_row, 1);
        Matrix Cy_temp(Cy_Nutations_full.n_row + Cy.n_row, 1);
        for (int m = 1; m <= Cx_Nutations_full.n_row; m++) {
            Cx_temp(m, 1) = Cx_Nutations_full(m, 1);
            Cy_temp(m, 1) = Cy_Nutations_full(m, 1);
        }
        for (int m = 1; m <= Cx.n_row; m++) {
            Cx_temp(Cx_Nutations_full.n_row + m, 1) = Cx(m, 1);
            Cy_temp(Cy_Nutations_full.n_row + m, 1) = Cy(m, 1);
        }
        Cx_Nutations_full = Cx_temp;
        Cy_Nutations_full = Cy_temp;
    }
    if (0 <= dt && dt <= 8) j = 0;
    else if (8 < dt && dt <= 16) j = 1;
    else if (16 < dt && dt <= 24) j = 2;
    else j = 3;
    Mjd0 = t1 + 8 * j;
    Matrix Cx_Nutations_sub(1, 10);
    Matrix Cy_Nutations_sub(1, 10);
    for (int k = 1; k <= 10; k++) {
        Cx_Nutations_sub(1, k) = Cx_Nutations_full(10 * j + k, 1);
        Cy_Nutations_sub(1, k) = Cy_Nutations_full(10 * j + k, 1);
    }
    Matrix Nutations = Cheb3D(Mjd_TDB, 10, Mjd0, Mjd0 + 8, Cx_Nutations_sub, Cy_Nutations_sub, zeros(1, 10));



	        std::cout << "AQUI LLEGASnutations";
			
			
    // Librations
    temp_aux = (929 - 899) / 10 + 1;
    temp = zeros(temp_aux, 1);
    for (int k = 1; k <= temp_aux; k++) {
        temp(k, 1) = (k - 1) * 10 + 899;
    }
    Matrix Cx_Librations(temp(2, 1) - temp(1, 1), 1);
    Matrix Cy_Librations(temp(3, 1) - temp(2, 1), 1);
    Matrix Cz_Librations(temp(4, 1) - temp(3, 1), 1);
    for (int col = temp(1, 1), idx = 1; col < temp(2, 1); col++, idx++) {
        Cx_Librations(idx, 1) = PCtemp(1, col);
    }
    for (int col = temp(2, 1), idx = 1; col < temp(3, 1); col++, idx++) {
        Cy_Librations(idx, 1) = PCtemp(1, col);
    }
    for (int col = temp(3, 1), idx = 1; col < temp(4, 1); col++, idx++) {
        Cz_Librations(idx, 1) = PCtemp(1, col);
    }
    Matrix Cx_Librations_full = Cx_Librations;
    Matrix Cy_Librations_full = Cy_Librations;
    Matrix Cz_Librations_full = Cz_Librations;
    for (int k = 1; k <= 3; k++) {
        temp = temp + 30;
        Matrix Cx(temp(2, 1) - temp(1, 1), 1);
        Matrix Cy(temp(3, 1) - temp(2, 1), 1);
        Matrix Cz(temp(4, 1) - temp(3, 1), 1);
        for (int col = temp(1, 1), idx = 1; col < temp(2, 1); col++, idx++) {
            Cx(idx, 1) = PCtemp(1, col);
        }
        for (int col = temp(2, 1), idx = 1; col < temp(3, 1); col++, idx++) {
            Cy(idx, 1) = PCtemp(1, col);
        }
        for (int col = temp(3, 1), idx = 1; col < temp(4, 1); col++, idx++) {
            Cz(idx, 1) = PCtemp(1, col);
        }
        Matrix Cx_temp(Cx_Librations_full.n_row + Cx.n_row, 1);
        Matrix Cy_temp(Cy_Librations_full.n_row + Cy.n_row, 1);
        Matrix Cz_temp(Cz_Librations_full.n_row + Cz.n_row, 1);
        for (int m = 1; m <= Cx_Librations_full.n_row; m++) {
            Cx_temp(m, 1) = Cx_Librations_full(m, 1);
            Cy_temp(m, 1) = Cy_Librations_full(m, 1);
            Cz_temp(m, 1) = Cz_Librations_full(m, 1);
        }
        for (int m = 1; m <= Cx.n_row; m++) {
            Cx_temp(Cx_Librations_full.n_row + m, 1) = Cx(m, 1);
            Cy_temp(Cy_Librations_full.n_row + m, 1) = Cy(m, 1);
            Cz_temp(Cz_Librations_full.n_row + m, 1) = Cz(m, 1);
        }
        Cx_Librations_full = Cx_temp;
        Cy_Librations_full = Cy_temp;
        Cz_Librations_full = Cz_temp;
    }
    if (0 <= dt && dt <= 8) j = 0;
    else if (8 < dt && dt <= 16) j = 1;
    else if (16 < dt && dt <= 24) j = 2;
    else j = 3;
    Mjd0 = t1 + 8 * j;
    Matrix Cx_Librations_sub(1, 10);
    Matrix Cy_Librations_sub(1, 10);
    Matrix Cz_Librations_sub(1, 10);
    for (int k = 1; k <= 10; k++) {
        Cx_Librations_sub(1, k) = Cx_Librations_full(10 * j + k, 1);
        Cy_Librations_sub(1, k) = Cy_Librations_full(10 * j + k, 1);
        Cz_Librations_sub(1, k) = Cz_Librations_full(10 * j + k, 1);
    }
    Matrix Librations = Cheb3D(Mjd_TDB, 10, Mjd0, Mjd0 + 8, Cx_Librations_sub, Cy_Librations_sub, Cz_Librations_sub);


	        std::cout << "AQUI LLEGASlibrations";

    // Position corrections
    const double EMRAT = 81.30056907419062; // DE430
    const double EMRAT1 = 1.0 / (1.0 + EMRAT);
    r_Earth = r_Earth - r_Moon * EMRAT1;
    r_Mercury = r_Mercury - r_Earth;
    r_Venus = r_Venus - r_Earth;
    r_Mars = r_Mars - r_Earth;
    r_Jupiter = r_Jupiter - r_Earth;
    r_Saturn = r_Saturn - r_Earth;
    r_Uranus = r_Uranus - r_Earth;
    r_Neptune = r_Neptune - r_Earth;
    r_Pluto = r_Pluto - r_Earth;
    r_Sun = r_Sun - r_Earth;
	
		        std::cout << "AQUI LLEGASFINALCPP";
				
				
				transpose(r_Mercury);
				transpose(r_Venus);
				transpose(r_Earth);
				transpose(r_Mars);
				transpose(r_Jupiter);
				transpose(r_Saturn);
				transpose(r_Uranus);
				transpose(r_Neptune);
				transpose(r_Pluto);
				transpose(r_Moon);
				transpose(r_Sun);
				
	
	
		
}


