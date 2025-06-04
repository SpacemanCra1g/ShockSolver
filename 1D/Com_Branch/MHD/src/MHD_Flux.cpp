#include "../include/DomainClass.hpp"

void Domain::Flux(double *Dest, double *P, int i, int destI) {
  double vx, vy, vz, d, p, Bx, By, Bz, TotP, E, Z;
  d = P[Tidx(DENSP, i)];
  vx = P[Tidx(VELX, i)];
  vy = P[Tidx(VELY, i)];
  vz = P[Tidx(VELZ, i)];
  p = P[Tidx(PRES, i)];
  Bx = P[Tidx(BX, i)];
  By = P[Tidx(BY, i)];
  Bz = P[Tidx(BZ, i)];
  TotP = p + .5 * (Bx * Bx + By * By + Bz * Bz) / d;

  E = p / ((GAMMA - 1) * d) + .5 * (vx * vx + vy * vy + vz * vz);
  Z = E + .5 * (Bx * Bx + By * By + Bz * Bz) / d;

  Dest[Tidx(DENS, destI)] = d * vx;
  Dest[Tidx(MOMX, destI)] = d * vx * vx + TotP - Bx * Bx;
  Dest[Tidx(MOMY, destI)] = d * vx * vy - Bx * By;
  Dest[Tidx(MOMZ, destI)] = d * vx * vz - Bx * Bz;
  Dest[Tidx(ENER, destI)] =
      vx * (d * Z + TotP) - Bx * (Bx * vx + vy * By + vz * Bz);
  Dest[Tidx(BX, destI)] = 0;
  Dest[Tidx(BY, destI)] = vx * By - Bx * vy;
  Dest[Tidx(BZ, destI)] = vx * Bz - Bx * vz;
}

void Domain::HLL_Flux(double *Dest, double *PrL, double *PrR, double SL,
                      double SR, int i) {

  double FluxL[NumVar], FluxR[NumVar];
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
