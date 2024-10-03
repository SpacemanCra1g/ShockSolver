#include "../include/VarConvert.hpp"

void Domain::Prims2Cons() {

  for (int i = 0; i < xDim; ++i) {
    MOMX[i] = DENS[i] * XVEL[i];
    ENERGY[i] = (XVEL[i] * XVEL[i] * DENS[i] / 2.0) + (PRES[i] / (GAMMA - 1.0));
  }
}

void Domain::Cons2Prim() {

  for (int i = 0; i < xDim; ++i) {
    XVEL[i] = XVEL[i] / DENS[i];
    PRES[i] = (GAMMA - 1.0) * (ENERGY[i] - XVEL[i] * XVEL[i] * DENS[i] / 2.0);
  }
}

void Domain::SolvePressure() {
  for (int i = 0; i < xDim; ++i) {
    PRES[i] = (GAMMA - 1.0) * (ENERGY[i] - XVEL[i] * XVEL[i] * DENS[i] / 2.0);
  }
}

void PrimConver(double *P, double *C) {

  C[Dens] = P[Dens];
  C[MomX] = P[Dens] * P[VelX];
  C[Ener] = (std::pow(P[VelX], 2) * P[Dens] / 2.0) + P[Pres] / (GAMMA - 1.0);
}

void ConConvert(double *C, double *P) {

  for (int i = 0; i < xDim; ++i) {
    P[Dens] = C[Dens];
    P[VelX] = C[MomX] / C[Dens];
    P[Pres] = (GAMMA - 1.0) * (C[Ener] - P[VelX] * P[VelX] * P[Dens] / 2.0);
  }
}

double HD_CS(double *P) { return std::sqrt(GAMMA * P[Pres] / P[Dens]); }

void SignalSpeed(double *P, double CS, double &CSL, double &CSR) {
  CSL = P[VelX] - CS;
  CSR = P[VelX] + CS;
}

void FillFlux(double P[], double F[]) {
  F[Dens] = P[Dens] * P[VelX];
  F[MomX] = P[Dens] * P[VelX] * P[VelX] + P[Pres];
  F[Ener] =
      (P[Pres] / (GAMMA - 1) + .5 * P[VelX] * P[VelX] * P[Dens] + P[Pres]) *
      P[VelX];
}
