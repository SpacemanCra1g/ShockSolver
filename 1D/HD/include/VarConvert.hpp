#ifndef VARCONVERT_H_
#define VARCONVERT_H_
#include "../include/DomainClass.hpp"

void PrimConver(double *P, double *C);

void ConConvert(double *C, double *P);

double HD_CS(double *P);

void SignalSpeed(double *P, double CS, double &CSL, double &CSR);

void FillFlux(double P[], double F[]);

#endif // VARCONVERT_H_
