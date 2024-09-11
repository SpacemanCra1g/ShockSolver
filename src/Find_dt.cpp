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

  for (int i = 0; i < xDim * yDim; ++i) {
    Cs[i] = std::sqrt(PRES[i] * GAMMA / DENS[i]);
  }
}

void Domain::Find_dt() {
  SolvePressure();
  Find_Cs();

  for (int i = 0; i < xDim * yDim; ++i) {
    Buffer[i] = dx / (std::fabs(XVEL[i]) + Cs[i]);
  }
  dt = FindMinimum(Buffer, xDim * yDim);

#if NDIMS > 1
  for (int i = 0; i < xDim * yDim; ++i) {
    Buffer[i] = dy / (std::fabs(YVEL[i] + Cs[i]));
  }

  dt = std::fmin(FindMinimum(Buffer, xDim * yDim), dt);
#endif

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
