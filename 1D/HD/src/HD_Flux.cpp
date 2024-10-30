#include "../include/DomainClass.hpp"

void Domain::HD_Flux(double *Dest, double *P, double *C, int i, int destI) {
  Dest[Tidx(DENS, destI)] = C[Tidx(MOMX, i)];
  Dest[Tidx(MOMX, destI)] =
      C[Tidx(MOMX, i)] * P[Tidx(VELX, i)] + P[Tidx(PRES, i)];

  Dest[Tidx(ENER, destI)] =
      (P[Tidx(PRES, i)] / (GAMMA - 1.0) +
       0.5 * P[Tidx(VELX, i)] * P[Tidx(VELX, i)] * P[Tidx(DENS, i)] +
       P[Tidx(PRES, i)]) *
      P[Tidx(VELX, i)];
}

void Domain::HD_HLL_Flux(double *Dest, double *PrL, double *CL, double *PrR,
                         double *CR, double SL, double SR, int i) {

  double FluxL[3], FluxR[3];

  FluxL[DENS] = CL[Tidx(MOMX, i)];
  FluxL[MOMX] = CL[Tidx(MOMX, i)] * PrL[Tidx(VELX, i)] + PrL[Tidx(PRES, i)];
  FluxL[ENER] =
      (PrL[Tidx(PRES, i)] / (GAMMA - 1.0) +
       0.5 * PrL[Tidx(VELX, i)] * PrL[Tidx(VELX, i)] * PrL[Tidx(DENS, i)] +
       PrL[Tidx(PRES, i)]) *
      PrL[Tidx(VELX, i)];

  FluxR[DENS] = CR[Tidx(MOMX, i + 1)];
  FluxR[MOMX] =
      CR[Tidx(MOMX, i + 1)] * PrR[Tidx(VELX, i + 1)] + PrR[Tidx(PRES, i + 1)];
  FluxR[ENER] = (PrR[Tidx(PRES, i + 1)] / (GAMMA - 1.0) +
                 0.5 * PrR[Tidx(VELX, i + 1)] * PrR[Tidx(VELX, i + 1)] *
                     PrR[Tidx(DENS, i + 1)] +
                 PrR[Tidx(PRES, i + 1)]) *
                PrR[Tidx(VELX, i + 1)];

  for (int var = 0; var < NumVar; ++var) {
    Dest[Tidx(var, i)] = (SR * FluxL[var] - SL * FluxR[var] +
                          SL * SR * (CR[Tidx(var, i + 1)] - CL[Tidx(var, i)])) /
                         (SR - SL);
  }
}
