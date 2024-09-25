#include "../include/DomainClass.hpp"

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
