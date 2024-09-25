#include "../include/DomainClass.hpp"

void Domain::Prims2Cons() {
  double g;
  double L;
  for (int i = 0; i < xDim; ++i) {
    L = Lorenz(i);
    g = DENSP[i] * Enthalpy(i) * L * L;

    MOMX[i] = L * XVEL[i];
    MOMY[i] = g * YVEL[i];
    MOMZ[i] = g * ZVEL[i];
    DENS[i] = g * DENSP[i];
    ENERGY[i] = g - PRES[i];
  }
}

void Domain::Cons2Prim() {
  double L;
  for (int i = 0; i < xDim; ++i) {
    Pressure(i);
    L = LorenzFromP(i);
    DENSP[i] = DENS[i] / L;
    XVEL[i] = MOMX[i] / (ENERGY[i] + PRES[i]);
    YVEL[i] = MOMY[i] / (ENERGY[i] + PRES[i]);
    ZVEL[i] = MOMZ[i] / (ENERGY[i] + PRES[i]);
  }
}

void Domain::SolvePressure() {
  for (int i = 0; i < xDim; ++i) {
    Pressure(i);
  }
}

double Domain::Lorenz(int x) {
  return std::pow(
      1.0 - XVEL[x] * XVEL[x] + YVEL[x] * YVEL[x] + ZVEL[x] * ZVEL[x], -0.5);
}

double Domain::LorenzFromP(int x) {
  double Ep = std::pow(ENERGY[x] + PRES[x], 2);
  return std::sqrt(
      Ep / (Ep + MOMX[x] * MOMX[x] + MOMY[x] * MOMY[x] + MOMZ[x] * MOMZ[x]));
}

double Domain::Enthalpy(int x) {
#if EOS == IdealGas
  return IDGas(x);
#endif
}

double Domain::EnthalpyFromCons(int x) {
#if EOS == IdealGas
  return IDGasFromCons(x);
#endif
}

double Domain::IDGasFromCons(int x) {
  return 1 + (GAMMA / (GAMMA - 1)) * (PRES[x] * Tau(x));
}

double Domain::IDGas(int x) {
  return 1 + (GAMMA / (GAMMA - 1)) * (PRES[x] / DENSP[x]);
}

double Domain::dh_dTau(int x) { return (GAMMA / (GAMMA - 1)) * PRES[x]; }

double Domain::dh_dP(int x) { return (GAMMA / (GAMMA - 1)) * Tau(x); }

double Domain::Tau(int x) { return (LorenzFromP(x) / DENS[x]); }

double Domain::F(int x) {
  return DENS[x] * EnthalpyFromCons(x) * LorenzFromP(x) - ENERGY[x] - PRES[x];
}

double Domain::dFp_dP(int x) {
  double L = LorenzFromP(x);
  double MNorm =
      std::pow(MOMX[x], 2) + std::pow(MOMY[x], 2) + std::pow(MOMZ[x], 2);
  return DENS[x] * L * dh_dP(x) -
         (MNorm * L * L * L / std::pow(ENERGY[x] + PRES[x], 3)) *
             (L * dh_dTau(x) + DENS[x] * EnthalpyFromCons(x)) -
         1.0;
}

double Domain::Newton(int x) { return PRES[x] - F(x) / dFp_dP(x); }

void Domain::Pressure(int x) {
  double PresLast;
  PRES[x] = 1.0;
  do {
    PresLast = PRES[x];
    PRES[x] = Newton(x);
    // std::cout << PRES[x] << std::endl;
  } while (std::fabs(PresLast - PRES[x]) > std::pow(10.0, -10.0));
}
