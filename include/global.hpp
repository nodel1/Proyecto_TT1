#ifndef _GLOBAL_
#define _GLOBAL_

#include "..\include\matrix.hpp"
#include <cmath>


extern Matrix eopdata;
extern Matrix Cnm;
extern Matrix Snm;
extern Matrix PC;

typedef struct {
    double Mjd_UTC, Mjd_TT;
    int n, m, sun, moon, planets;
}Param;

extern Param AuxParam;

void eop19620101(int c);
void GGM03S(int n);
void DE430Coeff(int f, int c);
void AuxParamInitialize();

#endif