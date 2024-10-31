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
    Cs[i] = std::sqrt(PRES[i] * GAMMA / DENS[i]);
  }
}

void Domain::Find_dt() {
  SolvePressure();
  Find_Cs();

  for (int i = 0; i < xDim; ++i) {
    if (count == 1) {
      // std::cout << Cs[i] << " Cell Number: " << i << std::endl;
    }
    Buffer[i] = dx / (std::fabs(XVEL[i]) + Cs[i]);
  }

  // if (count == 1) {
  //   exit(0);
  // }

  dt = FindMinimum(Buffer, xDim);
  dt *= CFL;

  if (T + dt > TN) {
    dt = TN - T;
  }

#if SlowStart
  if (dt > dt_sim) {
    dt = dt_sim;
    dt_sim *= 2.0;
  }
#endif
}
