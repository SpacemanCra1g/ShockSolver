#include "../include/DomainClass.hpp"

double FindMinimum(double *Array, const int size) {
  double Min = 1E100;
  for (int i = 0; i < size; ++i) {
    if (Min > Array[i]) {
      Min = Array[i];
    }
  }
  return Min;
}

void Domain::Find_Cs() {

  for (int i = 0; i < xDim; ++i) {
    // Cs[i] = HR_CS(i);
    Cs[i] = SRHD_CS(i);
  }
}

void Domain::Find_dt() {
  Cons2Prim();
  Find_Cs();
  double Lap;
  double nu;
  double Lor;

  for (int i = 0; i < xDim; ++i) {
    Lap = 1 - (XVEL[i] * XVEL[i] + YVEL[i] * YVEL[i] + ZVEL[i] * ZVEL[i]) *
                  Cs[i] * Cs[i];
    nu = 1 - XVEL[i] * XVEL[i] -
         Cs[i] * Cs[i] * (YVEL[i] * YVEL[i] + ZVEL[i] * ZVEL[i]);
    Lor = Lorenz(i);
    double Val = Lor * XVEL[i] * (1 - Cs[i] * Cs[i]);
    Buffer[i] =
        dx / std::max(std::fabs((Val - Cs[i] * std::sqrt(nu)) / (Lor * Lap)),
                      std::fabs((Val + Cs[i] * std::sqrt(nu)) / (Lor * Lap)));
  }
  dt = FindMinimum(Buffer, xDim);
  dt *= CFL;

  if (T + dt > TN) {
    dt = TN - T;
  }

#if SlowStart == True
  if (dt > dt_sim) {
    dt = dt_sim;
    dt_sim *= 2.0;
  }
#endif
}

double Domain::HD_CS(int i) { return std::sqrt(PRES[i] * GAMMA / DENS[i]); }

double Domain::SRHD_CS(int i) {
  return (std::pow(Tau(i), 2) / Enthalpy(i)) * dh_dTau(i) *
         (1.0 / (dh_dP(i) - Tau(i)));
}
