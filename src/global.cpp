#include "..\include\global.hpp"
#include "..\include\Mjday.hpp"
#include <stdio.h>
#include <string.h>
#include <cmath>

Matrix eopdata;
Matrix Cnm;
Matrix Snm;
Matrix PC;
Matrix obs; 

static const double Rad = 0.01745329251994329576923690768489;


void eop19620101(int c) {
 eopdata = zeros(13,c);
	
 FILE *fp = fopen("../DATA/eop19620101.txt","r");
	
 if (fp == NULL) {
  printf("Fail open eop19620101.txt file\n");
  exit(EXIT_FAILURE);
 }
	
 for (int j = 1; j <= c; j++) {
  fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &eopdata(1, j), 
  &eopdata(2, j),&eopdata(3, j),&eopdata(4, j),&eopdata(5, j),&eopdata(6, j),&eopdata(7, j),
  &eopdata(8, j),&eopdata(9, j),&eopdata(10, j),&eopdata(11, j),&eopdata(12, j),&eopdata(13, j));
 }
	        cout << "cargas eopdata";
 fclose(fp);
}



void GGM03S(int n){ 
    Cnm = zeros(n,n);
    Snm = zeros(n,n);
	std::cout << Cnm.n_row << " " << Cnm.n_column << " "  << Snm.n_row << " " << Snm.n_column << endl;
	
    FILE *fid = fopen("../DATA/GGM03S.txt","r");
	
    if(fid == NULL){
        printf("Fail open GGM03S.txt file\n");
        exit(EXIT_FAILURE);
    }

    double aux;
    for(int i=1; i<=n; i++){
        for(int j=1; j<=i; j++){
            fscanf(fid,"%lf %lf %lf %lf %lf %lf",&aux,&aux,&(Cnm(i,j)),&(Snm(i,j)),&aux,&aux);
        }
    }

    fclose(fid);
}

void DE430Coeff(int f, int c){
    PC = zeros(f,c);
	
    FILE *fid = fopen("../DATA/DE430Coeff.txt","r");
	
    if(fid == NULL){
        printf("Fail open DE430Coeff.txt file\n");
        exit(EXIT_FAILURE);
    }

    for(int i=1; i<=f; i++){
        for(int j=1; j<=c; j++){
            fscanf(fid,"%lf",&(PC(i,j)));
        }
    }

    fclose(fid);
}

Param AuxParam;

void AuxParamInitialize() {
    AuxParam.Mjd_UTC = 49746.1163541665;
    AuxParam.Mjd_TT = 49746.1170623147;
    AuxParam.n = 20;
    AuxParam.m = 20;
    AuxParam.sun = 1;
    AuxParam.moon = 1;
    AuxParam.planets = 1;
}



void readGEOS3(int nobs) {
    obs = zeros(nobs, 4);
    
    FILE *fid = fopen("../data/GEOS3.txt", "r");
    if (fid == NULL) {
        printf("Error al abrir GEOS3.txt\n");
        exit(EXIT_FAILURE);
    }

    char tline[256];
    char sub[11]; // Buffer para subcadenas (máximo 10 caracteres + '\0')
    int rows_read = 0;

    for (int i = 1; i <= nobs; i++) {
        if (!fgets(tline, sizeof(tline), fid) || strlen(tline) < 4) {
            break; // Fin de archivo o línea vacía
        }

        // Extraer año (posiciones 0-3)
        strncpy(sub, &tline[0], 4); sub[4] = '\0';
        int year = atoi(sub);

        // Extraer mes (posiciones 5-6)
        strncpy(sub, &tline[5], 2); sub[2] = '\0';
        int month = atoi(sub);

        // Extraer día (posiciones 8-9)
        strncpy(sub, &tline[8], 2); sub[2] = '\0';
        int day = atoi(sub);

        // Extraer hora (posiciones 11-12)
        strncpy(sub, &tline[11], 2); sub[2] = '\0';
        int hour = atoi(sub);

        // Extraer minuto (posiciones 14-15)
        strncpy(sub, &tline[14], 2); sub[2] = '\0';
        int min = atoi(sub);

        // Extraer segundo (posiciones 17-22, incluyendo decimales)
        strncpy(sub, &tline[17], 6); sub[6] = '\0';
        double sec = atof(sub);

        // Extraer azimut (posiciones 25-31)
        strncpy(sub, &tline[25], 7); sub[7] = '\0'; // Ajustado a 7 caracteres
        double az = atof(sub);

        // Extraer elevación (posiciones 34-40)
        strncpy(sub, &tline[34], 7); sub[7] = '\0'; // Ajustado a 7 caracteres
        double el = atof(sub);

        // Extraer distancia (posiciones 43-52)
        strncpy(sub, &tline[43], 9); sub[9] = '\0'; // Ajustado a 9 caracteres
        double dist = atof(sub);

        // Convertir a formato requerido
        obs(i, 1) = Mjday(year, month, day, hour, min, sec);
        obs(i, 2) = Rad * az; // Convertir a radianes
        obs(i, 3) = Rad * el; // Convertir a radianes
        obs(i, 4) = 1e3 * dist; // Convertir de km a m

        rows_read++;
    }

    fclose(fid);

    if (rows_read != nobs) {
        printf("Error: se esperaban %d observaciones, pero se leyeron %d.\n", nobs, rows_read);
        printf("Verifique que el archivo GEOS3.txt tenga el formato correcto y contenga suficientes datos.\n");
        exit(EXIT_FAILURE);
    }

    // Imprimir algunas observaciones para depuración
    printf("Datos cargados de GEOS3.txt:\n");
    for (int i = 1; i <= std::min(5, rows_read); i++) {
        printf("Obs %d: MJD = %.5f, Az = %.5f, El = %.5f, Dist = %.5f\n", 
               i, obs(i, 1), obs(i, 2), obs(i, 3), obs(i, 4));
    }
}





