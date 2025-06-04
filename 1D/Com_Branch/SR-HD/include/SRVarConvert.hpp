#ifndef GERVARCONVERT_H_
#define GERVARCONVERT_H_

#include "Parameters.h"
#include <cmath>

double IDGas(double *P);

double Enthalpy(double *P);

double Lorenz(double *P);

void PrimConvert(double *P, double Cons[5]);

double LorenzFromP(double *C, double PRES);

double dh_dTau(double PRES);

double Tau(double *C, double PRES);

double EnthalpyFromP(double *C, double PRES);

double dh_dP(double *C, double PRES);

double F(double *C, double PRES);

double dFp_dP(double *C, double PRES);

double Newton(double *C, double PRES);

double Pressure(double *C);

void ConConvert(double *C, double Prims[5]);

double SRH_CS(double *P);

void SignalSpeed(double *P, double CS, double &CSL, double &CSR);

#endif // GERVARCONVERT_H_
