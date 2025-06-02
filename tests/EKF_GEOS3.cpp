#include "..\include\matrix.hpp"
#include <cstdio>
#include <cmath>

#include "..\include\SAT_Const.hpp"
#include "..\include\global.hpp"

#include <iomanip>

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

static Matrix obs; // Declaración global estática
static const double Rad = 0.01745329251994329576923690768489; // Definición de Rad

Matrix AccelAdapter(double t, Matrix& Y) {
    return Accel(t, Y);  // Llama a Accel pasando Y (se hace una copia implícita)
}

void readGEOS3(int nobs) {
    obs = zeros(nobs, 4);
    
    FILE *fid = fopen("../data/GEOS3.txt", "r");
    if (fid == NULL) {
        printf("Error al abrir GEOS3.txt\n");
        exit(EXIT_FAILURE);
    }

    int year, month, day, hour, min;
    double sec, az, el, dist;
    
    for (int i = 1; i <= nobs; i++) {
        if (fscanf(fid, "%4d %2d %2d %2d %2d %6lf %8lf %7lf %10lf", 
                  &year, &month, &day, &hour, &min, &sec, 
                  &az, &el, &dist) != 9) {
            break;  // Error o fin de archivo
        }
        
        obs(i, 1) = Mjday(year, month, day, hour, min, sec);
        obs(i, 2) = Rad * az;
        obs(i, 3) = Rad * el;
        obs(i, 4) = 1e3 * dist;
    }

    fclose(fid);
}

Matrix& getObservations() {
    return obs;
}

int main() {
    // Cargar datos iniciales
    eop19620101(21413);  // Cargar parámetros de orientación terrestre
    GGM03S(181);         // Cargar coeficientes de gravedad (181x181)
    DE430Coeff(2285, 1020);  // Cargar coeficientes DE430
    AuxParamInitialize();  // Inicializar parámetros auxiliares
	
    int num_observations = 46;
    readGEOS3(num_observations);
	
    Matrix& observations = getObservations();
	
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

    Matrix r2(3, 1);
    Matrix v2(3, 1);
    r2(1, 1) = 6221397.62857869;
    r2(2, 1) = 2867713.77965738;
    r2(3, 1) = 3006155.98509949;
    v2(1, 1) = 4645.04725161806;
    v2(2, 1) = -2752.21591588204;
    v2(3, 1) = -7507.99940987031;

    Matrix Y0_apr(6, 1);
    Y0_apr(1, 1) = r2(1, 1);  // X
    Y0_apr(2, 1) = r2(2, 1);  // Y
    Y0_apr(3, 1) = r2(3, 1);  // Z
    Y0_apr(4, 1) = v2(1, 1);  // Vx
    Y0_apr(5, 1) = v2(2, 1);  // Vy
    Y0_apr(6, 1) = v2(3, 1);  // Vz

    // Época inicial
    double Mjd0 = Mjday(1995, 1, 29, 2, 38, 0);

    // Configurar parámetros para integración
    double Mjd_UTC = obs(9, 1);
    AuxParam.Mjd_UTC = Mjd_UTC;
    AuxParam.n = 20;
    AuxParam.m = 20;
    AuxParam.sun = true;
    AuxParam.moon = true;
    AuxParam.planets = true;
    int n_eqn = 6;

    // Integrar hacia atrás desde Mjd_UTC hasta Mjd0
    Matrix temp = DEInteg(AccelAdapter, 0, -(Mjd_UTC - Mjd0) * 86400.0, 1e-13, 1e-6, 6, Y0_apr);
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
        Matrix temp_variacional = DEInteg(AccelAdapter, 0, t - t_old, 1e-13, 1e-6, 42, yPhi);
        temporal2 = temp_variacional;
        yPhi = temporal2;

        // Extraer matriz de transición (sin assign_column, usando asignación manual)
        for (int j = 1; j <= 6; ++j) {
            for (int k = 1; k <= 6; ++k) {
                Phi(k, j) = yPhi(6 * j + k, 1);
            }
        }

        // Propagar estado
        Matrix temp_propagacion = DEInteg(AccelAdapter, 0, t - t_old, 1e-13, 1e-6, 6, Y_old);
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

        // Azimut y derivadas parciales
        AzElPa(s, Azim, Elev, dAds, dEds);
        dAdY = dAds * LT * U;
        dAdY = union_vector(dAdY, 1, 6); // Ajuste a union_vector

        // Actualización de medición (azimut)
        MeasUpdate(Y, obs(i, 2), Azim, sigma_az, dAdY, P, 6, K);

        // Elevación y derivadas parciales
        for (int j = 1; j <= 3; ++j) r(j, 1) = Y(j, 1);
        Ur = U * r;
        s = LT * (Ur - Rs);
        AzElPa(s, Azim, Elev, dAds, dEds);
        dEdY = dEds * LT * U;
        dEdY = union_vector(dEdY, 1, 6); // Ajuste a union_vector

        // Actualización de medición (elevación)
        MeasUpdate(Y, obs(i, 3), Elev, sigma_el, dEdY, P, 6, K);

        // Distancia y derivadas parciales
        for (int j = 1; j <= 3; ++j) r(j, 1) = Y(j, 1);
        Ur = U * r;
        s = LT * (Ur - Rs);
        double Dist = norm(s);
        dDds = transpose(s) / Dist;
        dDdY = dDds * LT * U;
        dDdY = union_vector(dDdY, 1, 6); // Ajuste a union_vector

        // Actualización de medición (distancia)
        MeasUpdate(Y, obs(i, 4), Dist, sigma_range, dDdY, P, 6, K);
    }

    // Integrar hacia atrás desde la última observación hasta la primera
    IERS(obs(46, 1), x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC);
    timediff(UT1_UTC, TAI_UTC, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);
    Mjd_TT = Mjd_UTC + TT_UTC / 86400.0;
    AuxParam.Mjd_UTC = Mjd_UTC;
    AuxParam.Mjd_TT = Mjd_TT;

    Matrix temp_Y0 = DEInteg(AccelAdapter, 0, -(obs(46, 1) - obs(1, 1)) * 86400.0, 1e-13, 1e-6, 6, Y);
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

    return 0;
}