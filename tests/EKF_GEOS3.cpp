#include "..\include\matrix.hpp"
#include <cstdio>
#include <cmath>

#include "..\include\SAT_Const.hpp"
#include "..\include\global.hpp"

#include <iomanip>
#include <chrono>

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
#include "..\include\NutMatrix.hpp"
#include "..\include\PrecMatrix.hpp"
#include "..\include\GMST.hpp"
#include "..\include\gast.hpp"
#include "..\include\MeasUpdate.hpp"
#include "..\include\G_AccelHarmonic.hpp"
#include "..\include\GHAMatrix.hpp"


#include "..\include\Accel.hpp"
#include "..\include\VarEqn.hpp"
#include "..\include\DEInteg.hpp"

#include <iomanip>


static const double Rad = 0.01745329251994329576923690768489; // Definición de Rad    lo mismo lo suyo era declararlo en satconst

Matrix AccelAdapter2(double t, Matrix& Y) {
Matrix result = Accel(t, Y);  // Llama a Accel y almacena el resultado
    return result;  // Devolver el resultado
}


Matrix VarEqnImpl(double t, Matrix& yPhi) {
    Matrix result = VarEqn(t, yPhi);  // Llama a la implementación interna y almacena el resultado
    return result;  // Devolver el resultado
}



Matrix& getObservations() {
    return obs;
}

int main() {
	
	    auto start = std::chrono::high_resolution_clock::now();
	

				
				
	cout << std::fixed << std::setprecision(15);
					
				
    // Cargar datos iniciales
    eop19620101(21413);  // Cargar parámetros de orientación terrestre
    GGM03S(181);         // Cargar coeficientes de gravedad (181x181)
    DE430Coeff(2285, 1020);  // Cargar coeficientes DE430
    AuxParamInitialize();  // Inicializar parámetros auxiliares


	
	
    int num_observations = 46;
    readGEOS3(num_observations);
	
    Matrix& observations = getObservations();
	
	// Verificar obs(9, 1)
    double obs_9_1 = obs(9, 1);
    printf("obs(9, 1) = %.5f\n", obs_9_1);
    if (obs_9_1 < 40000.0) { // Un MJD válido para 1995 debería ser alrededor de 49749
        printf("Error: obs(9, 1) tiene un valor no válido (%.5f). Se esperaba un MJD mayor a 40000.\n", obs_9_1);
        exit(EXIT_FAILURE);
    }

	
	
	
    // Parámetros de incertidumbre
    double sigma_range = 92.5;          
    double sigma_az = 0.0224 * Rad;    
    double sigma_el = 0.0139 * Rad;     

    // Estación Kaena Point
    double lat = Rad * 21.5748;         
    double lon = Rad * (-158.2706);     
    double alt = 300.20;                

    // Calcular posición de la estación
    Matrix pos = Position(lon, lat, alt);
    Matrix Rs = pos;

    // Seleccionar tiempos de observación
    double Mjd1 = obs(1, 1);
    double Mjd2 = obs(9, 1);
    double Mjd3 = obs(18, 1);

    Matrix r2= zeros(3, 1);
    Matrix v2= zeros(3, 1);
	Matrix aaaaaa= zeros(3, 1);
    r2(1, 1) = 6221397.62857869;
    r2(2, 1) = 2867713.77965738;
    r2(3, 1) = 3006155.98509949;
    v2(1, 1) = 4645.04725161806;
    v2(2, 1) = -2752.21591588204;
    v2(3, 1) = -7507.99940987031;

    Matrix Y0_apr = zeros(6, 1);
    Y0_apr(1, 1) = r2(1, 1);  // X
    Y0_apr(2, 1) = r2(2, 1);  // Y
    Y0_apr(3, 1) = r2(3, 1);  // Z
    Y0_apr(4, 1) = v2(1, 1);  // Vx
    Y0_apr(5, 1) = v2(2, 1);  // Vy
    Y0_apr(6, 1) = v2(3, 1);  // Vz

    // Época inicial
    double Mjd0 = Mjday(1995, 1, 29, 02, 38, 0);

    // Configurar parámetros para integración
    double Mjd_UTC = obs(9, 1);
    AuxParam.Mjd_UTC = Mjd_UTC;          
    AuxParam.n = 20;
    AuxParam.m = 20;
    AuxParam.sun = true;
    AuxParam.moon = true;
    AuxParam.planets = true;
    int n_eqn = 6;
	

				




// Imprimir el segundo argumento: t0 (double)
double t0 = 0;


// Imprimir el tercer argumento: tf (double)
double tf = -(Mjd_UTC - Mjd0) * 86400.0;                         //parece que esta aqui el error? en la resta del dei



// Imprimir el cuarto argumento: relerr (double)
double relerr = 1e-13;


// Imprimir el quinto argumento: abserr (double)
double abserr = 1e-6;





    // Integrar hacia atrás desde Mjd_UTC hasta Mjd0

    Matrix temp = DEInteg(AccelAdapter2, 0, tf, 1e-13, 1e-6, 6, Y0_apr); //parece que esta aqui el error? en la resta del dei del tf
    Matrix Y = temp;
	
	

	

	


    Matrix P = zeros(6, 6);
    for (int i = 1; i <= 3; ++i) {
        P(i, i) = 1e8;
    }
    for (int i = 4; i <= 6; ++i) {
        P(i, i) = 1e3;
    }

    // Calcular matriz de transformación local
    Matrix LT = LTC(lon, lat);

    // Variables para el bucle de medición
    Matrix yPhi = zeros(42, 1);  // Vector para integración variacional, 42x1
    Matrix temporal2 = zeros(42, 1); // Temporal
    Matrix Phi = zeros(6, 6);    // Matriz de transición de estado, 6x6

    // Tiempo
    double t = 0.0;              // Tiempo actual [s]
    double t_old;                // Tiempo anterior [s]

    // Parámetros de orientación terrestre y diferencias de tiempo
    double x_pole;               // Coordenada x del polo [rad]
    double y_pole;               // Coordenada y del polo [rad]
    double UT1_UTC;              // Diferencia UT1-UTC [s]
    double LOD;                  // Longitud del día [s]
    double dpsi;                 // Corrección de nutación en longitud [rad]
    double deps;                 // Corrección de nutación en oblicuidad [rad]
    double dx_pole;              // Corrección en x del polo [rad]
    double dy_pole;              // Corrección en y del polo [rad]
    double TAI_UTC;              // Diferencia TAI-UTC [s]
    double UT1_TAI;              // Diferencia UT1-TAI [s]
    double UTC_GPS;              // Diferencia UTC-GPS [s]
    double UT1_GPS;              // Diferencia UT1-GPS [s]
    double TT_UTC;               // Diferencia TT-UTC [s]
    double GPS_UTC;              // Diferencia GPS-UTC [s]
    double Mjd_TT;               // Fecha juliana modificada en TT
    double Mjd_UT1;              // Fecha juliana modificada en UT1
    double theta;                // Ángulo de rotación terrestre [rad]
    double Azim;                 // Azimut [rad]
    double Elev;                 // Elevación [rad]
	double Dist;

    // Matrices para el cálculo EKF
    Matrix Y_old(6, 1);          // Estado anterior (posición y velocidad), 6x1
    Matrix U(3, 3);              // Matriz de rotación terrestre, 3x3
    Matrix r(3, 1);              // Vector de posición, 3x1
    Matrix s(3, 1);              // Vector de posición topocéntrica, 3x1
    Matrix dAds(1, 3);           // Derivadas parciales de azimut respecto a s, 1x3
    Matrix dEds(1, 3);           // Derivadas parciales de elevación respecto a s, 1x3
    Matrix dAdY(1, 6);           // Derivadas parciales de azimut respecto a Y, 1x6
    Matrix K(6, 1);              // Ganancia de Kalman, 6x1
    Matrix dEdY(1, 6);           // Derivadas parciales de elevación respecto a Y, 1x6
    Matrix dDds(1, 3);           // Derivadas parciales de distancia respecto a s, 1x3
    Matrix dDdY(1, 6);           // Derivadas parciales de distancia respecto a Y, 1x6
    Matrix Ur(3, 1);             // Producto U * r

    // Bucle de medición (EKF)
    for (int i = 1; i <= 46; ++i) {
        t_old = t;
        Y_old = Y;

        // Incremento de tiempo y propagación
        Mjd_UTC = obs(i, 1);
        t = (Mjd_UTC - Mjd0) * 86400.0;

        IERS(Mjd_UTC, x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC);
        timediff(UT1_UTC, TAI_UTC, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);

        Mjd_TT = Mjd_UTC + TT_UTC / 86400.0;
        Mjd_UT1 = Mjd_TT + (UT1_UTC - TT_UTC) / 86400.0;
        AuxParam.Mjd_UTC = Mjd_UTC;
        AuxParam.Mjd_TT = Mjd_TT;

        // Construir vector yPhi (estado + matriz de transición)
        for (int ii = 1; ii <= 6; ++ii) {
            yPhi(ii, 1) = Y_old(ii, 1);
            for (int j = 1; j <= 6; ++j) {
                yPhi(6 * j + ii, 1) = (ii == j) ? 1.0 : 0.0;
            }
        }

        // Integrar ecuaciones variacionales

        Matrix temp_variacional = DEInteg(&VarEqnImpl, 0, t - t_old, 1e-13, 1e-6, 42, yPhi);
        temporal2 = temp_variacional;
        yPhi = temporal2;


        // Extraer matriz de transición (sin assign_column, usando asignación manual)
        for (int j = 1; j <= 6; ++j) {
            for (int k = 1; k <= 6; ++k) {
                Phi(k, j) = yPhi(6 * j + k, 1);
            }
        }

        // Propagar estado

        Matrix temp_propagacion = DEInteg(AccelAdapter2, 0, t - t_old, 1e-13, 1e-6, 6, Y_old);
        Y = temp_propagacion;

				

        // Coordenadas topocéntricas
        theta = GMST(Mjd_UT1);
        Matrix temp_U = R_z(theta);
		

						
						
        U = temp_U;
        for (int j = 1; j <= 3; ++j) r(j, 1) = Y(j, 1);
        Ur = U * r;
        s = LT * (Ur - Rs);
		


        // Actualización de tiempo (propagación de covarianza)
        Matrix temp_P = TimeUpdate(P, Phi);
        P = temp_P;
		

								

AzElPa(s, Azim, Elev, dAds, dEds);




dAdY = dAds * LT * U;


// Extend dAdY to 1x6 with zero velocity derivatives
Matrix full_dAdY(1, 6);
for (int j = 1; j <= 3; j++) {
    full_dAdY(1, j) = dAdY(1, j); // Copy position derivatives
}
for (int j = 4; j <= 6; j++) {
    full_dAdY(1, j) = 0.0; // Velocity derivatives are zero
}
dAdY = full_dAdY;


// Remove union_vector since we manually reshaped to 1x6




// Actualización de medición (azimut)
MeasUpdate(Y, obs(i, 2), Azim, sigma_az, dAdY, P, 6, K);


        // Elevación y derivadas parciales
for (int j = 1; j <= 3; ++j) r(j, 1) = Y(j, 1);
Ur = U * r;
s = LT * (Ur - Rs);
AzElPa(s, Azim, Elev, dAds, dEds);
dEdY = dEds * LT * U;


// Extend dEdY to 1x6 with zero velocity derivatives
Matrix full_dEdY(1, 6);
for (int j = 1; j <= 3; j++) {
    full_dEdY(1, j) = dEdY(1, j); // Copy position derivatives
}
for (int j = 4; j <= 6; j++) {
    full_dEdY(1, j) = 0.0; // Velocity derivatives are zero
}
dEdY = full_dEdY;






// Actualización de medición (elevación)
MeasUpdate(Y, obs(i, 3), Elev, sigma_el, dEdY, P, 6, K);


// Distancia y derivadas parciales
for (int j = 1; j <= 3; ++j) r(j, 1) = Y(j, 1);
Ur = U * r;
s = LT * (Ur - Rs);
Dist = norm(s);
dDds = transpose(s) / Dist;
dDdY = dDds * LT * U;


// Extend dDdY to 1x6 with zero velocity derivatives
Matrix full_dDdY(1, 6);
for (int j = 1; j <= 3; j++) {
    full_dDdY(1, j) = dDdY(1, j); // Copy position derivatives
}
for (int j = 4; j <= 6; j++) {
    full_dDdY(1, j) = 0.0; // Velocity derivatives are zero
}
dDdY = full_dDdY;





// Actualización de medición (distancia)
MeasUpdate(Y, obs(i, 4), Dist, sigma_range, dDdY, P, 6, K);

    }



    // Integrar hacia atrás desde la última observación hasta la primera
    IERS(obs(46, 1), x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC);
    timediff(UT1_UTC, TAI_UTC, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);
    Mjd_TT = Mjd_UTC + TT_UTC / 86400.0;
    AuxParam.Mjd_UTC = Mjd_UTC;
    AuxParam.Mjd_TT = Mjd_TT;



    Matrix temp_Y0 = DEInteg(AccelAdapter2, 0, -(obs(46, 1) - obs(1, 1)) * 86400.0, 1e-13, 1e-6, 6, Y);
    Matrix Y0 = temp_Y0;




    // Estado verdadero para comparación
    Matrix Y_true(6, 1);
    Y_true(1, 1) = 5753.173e3;  // Posición X [m]
    Y_true(2, 1) = 2673.361e3;  // Posición Y [m]
    Y_true(3, 1) = 3440.304e3;  // Posición Z [m]
    Y_true(4, 1) = 4.324207e3;  // Velocidad X [m/s]
    Y_true(5, 1) = -1.924299e3; // Velocidad Y [m/s]
    Y_true(6, 1) = -5.728216e3; // Velocidad Z [m/s]

    // Imprimir errores
    std::cout << "\nError de Estimación de Posición\n";
    std::cout << std::fixed << std::setprecision(10);
    std::cout << "dX " << std::setw(10) << Y0(1, 1) - Y_true(1, 1) << " [m]\n";
    std::cout << "dY " << std::setw(10) << Y0(2, 1) - Y_true(2, 1) << " [m]\n";
    std::cout << "dZ " << std::setw(10) << Y0(3, 1) - Y_true(3, 1) << " [m]\n";
    std::cout << "\nError de Estimación de Velocidad\n";
    std::cout << "dVx " << std::setw(8) << Y0(4, 1) - Y_true(4, 1) << " [m/s]\n";
    std::cout << "dVy " << std::setw(8) << Y0(5, 1) - Y_true(5, 1) << " [m/s]\n";
    std::cout << "dVz " << std::setw(8) << Y0(6, 1) - Y_true(6, 1) << " [m/s]\n";
	
	
	auto end = std::chrono::high_resolution_clock::now(); // Fin

    std::chrono::duration<double> duration = end - start;
    std::cout << "Tiempo de ejecución (con impresiones por pantalla, revisar): " << duration.count() << " segundos" << std::endl;
	
	

    return 0;
}