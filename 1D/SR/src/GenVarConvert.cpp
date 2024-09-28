#include "../include/GenVarConvert.hpp"

double IDGas(double *P) {
  return 1 + (GAMMA / (GAMMA - 1.0)) * (P[Pres] / P[DensP]);
}

double Enthalpy(double *P) {
#if EOS == IdealGas
  return IDGas(P);
#endif
}

double Lorenz(double *P) {
  return std::pow(
      1.0 - (P[VelX] * P[VelX] + P[VelY] * P[VelY] + P[VelZ] * P[VelZ]), -0.5);
}

void Prims2Cons(double *P, double Cons[5]) {
  double g;
  double L;

  L = Lorenz(P);
  g = P[DensP] * Enthalpy(P) * L * L;

  Cons[0] = L * P[0];

  for (int i = 1; i < 4; ++i) {
    Cons[i] = g * P[i];
  }
  Cons[4] = g - P[4];
}

double LorenzFromP(double *C, double PRES) {
  double Ep = std::pow(C[Ener] + PRES, 2);
  double Norm = 0.0;

  for (int i = MomX; i <= MomZ; ++i) {
    Norm += C[i] * C[i];
  }

  return 1.0 / std::sqrt(1 - Norm / Ep);
}

double dh_dTau(double PRES) { return (GAMMA / (GAMMA - 1)) * PRES; }

double Tau(double *C, double PRES) { return (LorenzFromP(C, PRES) / C[Dens]); }

double EnthalpyFromP(double *C, double PRES) {
#if EOS == IdealGas
  return 1 + (GAMMA / (GAMMA - 1.0)) * (PRES * Tau(C, PRES));
#endif
}

double dh_dP(double *C, double PRES) {
  return (GAMMA / (GAMMA - 1)) * Tau(C, PRES);
}

double F(double *C, double PRES) {
  return C[Dens] * EnthalpyFromP(C, PRES) * LorenzFromP(C, PRES) - C[Ener] -
         PRES;
}

double dFp_dP(double *C, double PRES) {
  double L = LorenzFromP(C, PRES);
  double MNorm =
      std::pow(C[MomX], 2) + std::pow(C[MomY], 2) + std::pow(C[MomZ], 2);
  return C[Dens] * L * dh_dP(C, PRES) -
         (MNorm * L * L * L / std::pow(C[Ener] + PRES, 3)) *
             (L * dh_dTau(PRES) + C[Dens] * EnthalpyFromP(C, PRES)) -
         1.0;
}

double Newton(double *C, double PRES) {
  return PRES - F(C, PRES) / dFp_dP(C, PRES);
}

double Pressure(double *C) {
  double PresLast;
  double PRES;
  PRES = 1.0;
  do {
    PresLast = PRES;
    PRES = Newton(C, PRES);
  } while (std::fabs(PresLast - PRES) > std::pow(10.0, -12.0));

  return PRES;
}

void Cons2Prim(double *C, double Prims[5]) {
  {
    double L;
    for (int i = 0; i < xDim; ++i) {
      Prims[Pres] = Pressure(C);
      L = LorenzFromP(C, Prims[Pres]);
      Prims[DensP] = C[Dens] / L;

      for (int i = VelX; i <= VelZ; ++i) {
        Prims[i] = C[i] / (C[Ener] + Prims[Pres]);
      }
    }
  }
}

double SRHD_CS(double *C, double *P) {
  return (std::pow(Tau(C, P[Pres]), 2) / Enthalpy(P)) * dh_dTau(P[Pres]) *
         (1.0 / (dh_dP(C, P[Pres]) - Tau(C, P[Pres])));
}

void SignalSpeed(double *P, double CS, double &CSL, double &CSR) {
  double Norm = 0.0;
  for (int i = VelX; i <= VelZ; ++i) {
    Norm += P[i] * P[i];
  }
  double Lap = 1.0 - Norm * CS;
  double nu =
      1.0 - P[VelX] * P[VelX] - CS * (P[VelY] * P[VelY] + P[VelZ] * P[VelZ]);
  double Lor = Lorenz(P);
  double Val = Lor * P[VelX] * (1 - CS);

  CSL = (Val - std::sqrt(CS * nu)) / (Lor * Lap);
  CSR = (Val + std::sqrt(nu * CS)) / (Lor * Lap);
}
