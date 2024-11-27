#include "../include/DomainClass.hpp"
#include <iostream>

double IDGas(double *P) {
  return 1 + (GAMMA / (GAMMA - 1.0)) * (P[PRES] / P[DENSP]);
}

double Enthalpy(double *P) {
#if EOS == IDEAL
  return IDGas(P);
#endif
}

double Lorenz(double *P) {

  double Norm = 0.0;
  for (int var = VELX; var <= VELZ; var++) {
    Norm += std::pow(P[var], 2);
    if (Norm >= 1.0) {
      std::cout << "Vel excedes 1 " << Norm << std::endl;
      exit(0);
    }
  }

  return std::pow(1.0 - Norm, -0.5);
}

void PrimConvert(double *P, double *Cons) {
  double g;
  double L;

  L = Lorenz(P);
  g = P[DENSP] * Enthalpy(P) * L * L;

  Cons[DENS] = L * P[DENSP];

  for (int i = VELX; i <= VELZ; ++i) {
    Cons[i] = g * P[i];
  }
  Cons[4] = g - P[PRES];
}

double LorenzFromP(double *C, double Pres) {
  double Ep = std::pow(C[ENER] + Pres, 2);
  double Norm = 0.0;

  for (int i = MOMX; i <= MOMZ; ++i) {
    Norm += C[i] * C[i];
  }
  // if (Norm >= 1.0) {
  // std::cout << "Vel excedes 1 " << Norm << std::endl;
  // exit(0);
  // }

  return 1.0 / std::sqrt(1 - Norm / Ep);
}

double dh_dTau(double Pres) { return (GAMMA / (GAMMA - 1)) * Pres; }

double Tau(double *C, double Pres) { return (LorenzFromP(C, Pres) / C[DENS]); }

double EnthalpyFromP(double *C, double Pres) {
#if EOS == IDEAL
  return 1 + (GAMMA / (GAMMA - 1.0)) * (Pres * Tau(C, Pres));
#endif
}

double dh_dP(double *C, double Pres) {
  return (GAMMA / (GAMMA - 1)) * Tau(C, Pres);
}

double F(double *C, double Pres) {
  return C[DENS] * EnthalpyFromP(C, Pres) * LorenzFromP(C, Pres) - C[ENER] -
         Pres;
}

double dFp_dP(double *C, double Pres) {
  double L = LorenzFromP(C, Pres);
  double MNorm = 0.0;

  for (int var = MOMX; var <= MOMZ; ++var) {
    MNorm += std::pow(C[var], 2);
  }

  return C[DENS] * L * dh_dP(C, Pres) -
         (MNorm * L * L * L / std::pow(C[ENER] + Pres, 3)) *
             (L * dh_dTau(Pres) + C[DENS] * EnthalpyFromP(C, Pres)) -
         1.0;
}

double Newton(double *C, double Pres) {
  double var = Pres - F(C, Pres) / dFp_dP(C, Pres);
  if (var <= 0.0) {
    var = std::pow(10.0, -10);
  }
  return var;
}

double Pressure(double *C) {
  double PresLast;
  double Pres;
  int count = 0;
  double check = -1.0;
  Pres = pow(10.0, check);
  do {
    PresLast = Pres;
    Pres = Newton(C, Pres);
    count += 1;
    if (count == 30) {
      std::cout << "Failure in Newton's method Count Exceded" << std::endl;
      return -1.0;
    }

  } while (std::fabs(PresLast - Pres) > std::pow(10.0, -9.0));

  return Pres;
}

int Domain::NaiveNewton(double *Uin, double *Uout, int i) {

  double L;
  double C[5], P[5];
  for (int var = 0; var < NumVar; ++var) {
    C[var] = Uin[Tidx(var, i)];
  }

  P[PRES] = Pressure(C);
  if (P[PRES] < 0.0) {
    return 1;
  }
  L = LorenzFromP(C, P[PRES]);
  P[DENSP] = C[DENS] / L;

  for (int i = VELX; i <= VELZ; ++i) {
    P[i] = C[i] / (C[ENER] + P[PRES]);
  }
  for (int var = 0; var < NumVar; ++var) {
    Uout[Tidx(var, i)] = P[var];
  }
  return 0;
}

// double SRH_CS(double *P) { return GAMMA * P[Pres] / (P[DensP] * Enthalpy(P));
// }

// void SignalSpeed(double *P, double CS, double &CSL, double &CSR) {
//   double Norm = 0.0;
//   for (int i = VelX; i <= VelZ; ++i) {
//     Norm += P[i] * P[i];
//   }

//   double sroot = std::sqrt(
//       CS * (1 - P[VelX] * P[VelX] -
//             (P[VelY] * P[VelY] + P[VelZ] * P[VelZ]) * CS * (1 - Norm)));

//   CSR = (P[VelX] * (1 - CS) + sroot) / (1 - Norm * CS);
//   CSL = (P[VelX] * (1 - CS) - sroot) / (1 - Norm * CS);
// }
