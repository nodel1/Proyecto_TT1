#include "..\include\global.hpp"
#include "..\include\Mjday.hpp"
#include <stdio.h>
#include <string.h>
#include <cmath>
#include <iomanip>



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
    std::cout << std::fixed << std::setprecision(15);  // Para otras salidas con cout
    obs = zeros(nobs, 4);
    
    FILE *fid = fopen("../data/GEOS3.txt", "r");
    if (fid == NULL) {
        printf("Error al abrir GEOS3.txt\n");
        exit(EXIT_FAILURE);
    }

    char tline[256];
    int rows_read = 0;

    for (int i = 1; i <= nobs; i++) {
        if (!fgets(tline, sizeof(tline), fid) || strlen(tline) < 4) {
            break;
        }

        int year, month, day, hour, min;
        double sec, az, el, dist;

        // Parsear la línea con sscanf para el formato YYYY/MM/DD HH:MM:SS.SS
        int fields = sscanf(tline, "%d/%d/%d %d:%d:%lf %lf %lf %lf",
                            &year, &month, &day, &hour, &min, &sec, &az, &el, &dist);

        if (fields != 9) {
            printf("Error al parsear la fila %d: se leyeron %d campos, se esperaban 9.\n", i, fields);
            printf("Linea: %s", tline);
            fclose(fid);
            exit(EXIT_FAILURE);
        }

        // Debug para la fila 9
        if (i == 9) {
            printf("Fila 9 antes de Mjday:\n");
            printf("  Linea leida: %s", tline);
            printf("  Year: %d\n", year);
            printf("  Month: %d\n", month);
            printf("  Day: %d\n", day);
            printf("  Hour: %d\n", hour);
            printf("  Min: %d\n", min);
            printf("  Sec: %.10f\n", sec);
            printf("  Az: %.10f\n", az);
            printf("  El: %.10f\n", el);
            printf("  Dist: %.10f\n", dist);
        }

        // Guardar en matriz
        obs(i, 1) = Mjday(year, month, day, hour, min, sec);
        obs(i, 2) = Rad * az;
        obs(i, 3) = Rad * el;
        obs(i, 4) = 1e3 * dist;

        rows_read++;
    }

    fclose(fid);

    if (rows_read != nobs) {
        printf("Error: se esperaban %d observaciones, pero se leyeron %d.\n", nobs, rows_read);
        printf("Verifique que el archivo GEOS3.txt tenga el formato correcto y contenga suficientes datos.\n");
        exit(EXIT_FAILURE);
    }

    // Mostrar con 10 decimales todas las columnas
    printf("Datos cargados de GEOS3.txt:\n");
    for (int i = 1; i <= std::min(5, rows_read); i++) {
        printf("Obs %2d: MJD = %.10f, Az = %.10f, El = %.10f, Dist = %.10f\n", 
               i, obs(i, 1), obs(i, 2), obs(i, 3), obs(i, 4));
    }
}





