#include "../include/DomainClass.hpp"

void Domain::SR_Flux(double *Dest, double *P, double *C, int i, int destI) {
  Dest[Tidx(DENS, destI)] = C[Tidx(DENS, i)] * P[Tidx(VELX, i)];
  Dest[Tidx(MOMX, destI)] =
      C[Tidx(MOMX, i)] * P[Tidx(VELX, i)] + P[Tidx(PRES, i)];
  Dest[Tidx(MOMY, destI)] = C[Tidx(MOMY, i)] * P[Tidx(VELX, i)];
  Dest[Tidx(MOMZ, destI)] = C[Tidx(MOMZ, i)] * P[Tidx(VELX, i)];
  Dest[Tidx(ENER, destI)] = C[Tidx(MOMX, i)];
}

void Domain::SR_HLL_Flux(double *Dest, double *PrL, double *CL, double *PrR,
                         double *CR, double SL, double SR, int i) {

  double FluxL[5], FluxR[5];

  FluxL[DENS] = CL[Tidx(DENS, i)] * PrL[Tidx(VELX, i)];
  FluxL[MOMX] = CL[Tidx(MOMX, i)] * PrL[Tidx(VELX, i)] + PrL[Tidx(PRES, i)];
  FluxL[MOMY] = CL[Tidx(MOMY, i)] * PrL[Tidx(VELX, i)];
  FluxL[MOMZ] = CL[Tidx(MOMZ, i)] * PrL[Tidx(VELX, i)];
  FluxL[ENER] = CL[Tidx(MOMX, i)];

  FluxR[DENS] = CR[Tidx(DENS, i + 1)] * PrR[Tidx(VELX, i + 1)];
  FluxR[MOMX] =
      CR[Tidx(MOMX, i + 1)] * PrR[Tidx(VELX, i + 1)] + PrR[Tidx(PRES, i + 1)];
  FluxR[MOMY] = CR[Tidx(MOMY, i + 1)] * PrR[Tidx(VELX, i + 1)];
  FluxR[MOMZ] = CR[Tidx(MOMZ, i + 1)] * PrR[Tidx(VELX, i + 1)];
  FluxR[ENER] = CR[Tidx(MOMX, i + 1)];

  for (int var = 0; var < NumVar; ++var) {
    Dest[Tidx(var, i)] = (SR * FluxL[var] - SL * FluxR[var] +
                          SL * SR * (CR[Tidx(var, i + 1)] - CL[Tidx(var, i)])) /
                         (SR - SL);
  }
}
