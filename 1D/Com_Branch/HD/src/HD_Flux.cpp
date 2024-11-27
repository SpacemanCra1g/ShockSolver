#include "../include/DomainClass.hpp"

void Domain::Flux(double *Dest, double *P, int i, int destI) {
  double vx, d, p;
  d = P[Tidx(DENSP, i)];
  vx = P[Tidx(VELX, i)];
  p = P[Tidx(PRES, i)];

  Dest[Tidx(DENS, destI)] = d * vx;
  Dest[Tidx(MOMX, destI)] = d * vx * vx + p;
  Dest[Tidx(ENER, destI)] =
      vx * (d * (0.5 * vx * vx + p / ((GAMMA - 1.0) * d)) + p);
}

void Domain::HLL_Flux(double *Dest, double *PrL, double *PrR, double SL,
                      double SR, int i) {

  double FluxL[3], FluxR[3];
  double dL, vxL, pL;
  double dR, vxR, pR, scalar, EL, ER;

  dL = PrL[Tidx(DENSP, i)];
  vxL = PrL[Tidx(VELX, i)];
  pL = PrL[Tidx(PRES, i)];

  dR = PrR[Tidx(DENSP, i + 1)];
  vxR = PrR[Tidx(VELX, i + 1)];
  pR = PrR[Tidx(PRES, i + 1)];

  FluxL[DENS] = dL * vxL;
  FluxL[MOMX] = dL * vxL * vxL + pL;
  EL = dL * (0.5 * vxL * vxL + pL / ((GAMMA - 1.0) * dL));
  FluxL[ENER] = vxL * (EL + pL);

  FluxR[DENS] = dR * vxR;
  FluxR[MOMX] = dR * vxR * vxR + pR;
  // std::cout << dR << std::endl;
  ER = dR * (0.5 * vxR * vxR + pR / ((GAMMA - 1.0) * dR));
  FluxR[ENER] = vxR * (ER + pR);

  for (int var = 0; var < NumVar; ++var) {
    Dest[Tidx(var, i)] = (SR * FluxL[var] - SL * FluxR[var]) / (SR - SL);
  }
  scalar = SL * SR / (SR - SL);
  Dest[Tidx(DENS, i)] += scalar * (dR - dL);
  Dest[Tidx(MOMX, i)] += scalar * (dR * vxR - dL * vxL);
  Dest[Tidx(ENER, i)] += scalar * (ER - EL);
}
