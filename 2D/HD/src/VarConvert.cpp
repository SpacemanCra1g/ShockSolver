#include "../include/DomainClass.hpp"

void Domain::Prims2Cons() {
  // #pragma omp simd
  //   {
  for (int i = 0; i < xDim * yDim; ++i) {
    MOMX[i] = DENS[i] * XVEL[i];

#if NDIMS > 1
    // Energy and pres conversions are wrong in 2D
    ENERGY[i] = PRES[i] / (GAMMA - 1) +
                .5 * DENS[i] * ((XVEL[i] * XVEL[i]) + YVEL[i] * YVEL[i]);
    MOMY[i] = DENS[i] * YVEL[i];
#else
    ENERGY[i] = (XVEL[i] * XVEL[i] * DENS[i] / 2.0) + (PRES[i] / (GAMMA - 1.0));

#endif
    //
  }
}

void Domain::Cons2Prim() {
  // #pragma omp simd
  //   {
  for (int i = 0; i < xDim * yDim; ++i) {
    XVEL[i] = XVEL[i] / DENS[i];

#if NDIMS > 1
    // See above
    MOMY[i] = YVEL[i] / DENS[i];
    PRES[i] =
        (GAMMA - 1) *
        (ENERGY[i] - 0.5 * (MOMX[i] * MOMX[i] + MOMY[i] * MOMY[i]) / DENS[i]);
#else

    PRES[i] = (GAMMA - 1.0) * (ENERGY[i] - XVEL[i] * XVEL[i] * DENS[i] / 2.0);
#endif
  }
  // }
}

void Domain::SolvePressure() {
  for (int i = 0; i < xDim * yDim; ++i) {
#if NDIMS > 1
    PRES[i] =
        (GAMMA - 1) *
        (ENERGY[i] - 0.5 * (MOMX[i] * MOMX[i] + MOMY[i] * MOMY[i]) / DENS[i]);
#else
    PRES[i] = (GAMMA - 1.0) * (ENERGY[i] - XVEL[i] * XVEL[i] * DENS[i] / 2.0);
#endif
  }
}

void Cell::Cons2Prims() {
  *XVEL = *XVEL / *DENS;
  *MOMY = *YVEL / *DENS;
  *PRES =

      (GAMMA - 1) * (*ENERGY - 0.5 * (*MOMX * *MOMX + *MOMY * *MOMY) / *DENS);
}

void Cell::Prims2Cons() {
  *MOMX = *DENS * *XVEL;
  *MOMY = *DENS * *YVEL;
  *ENERGY =

      *PRES / (GAMMA - 1) + .5 * *DENS * ((*XVEL * *XVEL) + *YVEL * *YVEL);
}

double Cell::GetPres() {
  *PRES =

      (GAMMA - 1) * (*ENERGY - 0.5 * (*MOMX * *MOMX + *MOMY * *MOMY) / *DENS);

  return *PRES;
}
