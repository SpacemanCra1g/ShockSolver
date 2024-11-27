#include "../include/DomainClass.hpp"

void Domain::Flux(double *Dest, double *P, int i, int destI) {
  double vx, vy, vz, d, p, lor, val, h;
  d = P[Tidx(DENSP, i)];
  vx = P[Tidx(VELX, i)];
  vy = P[Tidx(VELY, i)];
  vz = P[Tidx(VELZ, i)];
  p = P[Tidx(PRES, i)];
  lor = std::pow(1.0 - (vx * vx + vy * vy + vz * vz), -0.5);
#if EOS == IDEAL
  h = 1.0 + (GAMMA / (GAMMA - 1.0)) * p / d;
#endif
  val = d * h * lor * lor;

  Dest[Tidx(DENS, destI)] = lor * d * vx;
  Dest[Tidx(MOMX, destI)] = val * vx * vx + p;
  Dest[Tidx(MOMY, destI)] = val * vy * vx;
  Dest[Tidx(MOMZ, destI)] = val * vz * vx;
  Dest[Tidx(ENER, destI)] = val * vx;
}

void Domain::HLL_Flux(double *Dest, double *PrL, double *PrR, double SL,
                      double SR, int i) {

  double FluxL[5], FluxR[5];
  double dL, vxL, vyL, vzL, pL;
  double dR, vxR, vyR, vzR, pR;
  double lorL, hL, valL, lorR, hR, valR, scalar;

  dL = PrL[Tidx(DENSP, i)];
  vxL = PrL[Tidx(VELX, i)];
  vyL = PrL[Tidx(VELY, i)];
  vzL = PrL[Tidx(VELZ, i)];
  pL = PrL[Tidx(PRES, i)];

  dR = PrR[Tidx(DENSP, i + 1)];
  vxR = PrR[Tidx(VELX, i + 1)];
  vyR = PrR[Tidx(VELY, i + 1)];
  vzR = PrR[Tidx(VELZ, i + 1)];
  pR = PrR[Tidx(PRES, i + 1)];

  lorL = std::pow(1.0 - (vxL * vxL + vyL * vyL + vzL * vzL), -0.5);
  lorR = std::pow(1.0 - (vxR * vxR + vyR * vyR + vzR * vzR), -0.5);
#if EOS == IDEAL
  hL = 1.0 + (GAMMA / (GAMMA - 1.0)) * pL / dL;
  hR = 1.0 + (GAMMA / (GAMMA - 1.0)) * pR / dR;
#endif
  valL = dL * hL * lorL * lorL;
  valR = dR * hR * lorR * lorR;

  FluxL[DENS] = lorL * dL * vxL;
  FluxL[MOMX] = valL * vxL * vxL + pL;
  FluxL[MOMY] = valL * vyL * vxL;
  FluxL[MOMZ] = valL * vzL * vxL;
  FluxL[ENER] = valL * vxL;

  FluxR[DENS] = lorR * dR * vxR;
  FluxR[MOMX] = valR * vxR * vxR + pR;
  FluxR[MOMY] = valR * vyR * vxR;
  FluxR[MOMZ] = valR * vzR * vxR;
  FluxR[ENER] = valR * vxR;

  for (int var = 0; var < NumVar; ++var) {
    Dest[Tidx(var, i)] = (SR * FluxL[var] - SL * FluxR[var]) / (SR - SL);
  }
  scalar = SL * SR / (SR - SL);
  Dest[Tidx(DENS, i)] += scalar * (lorR * dR - lorL * dL);
  Dest[Tidx(MOMX, i)] += scalar * (valR * vxR - valL * vxL);
  Dest[Tidx(MOMY, i)] += scalar * (valR * vyR - valL * vyL);
  Dest[Tidx(MOMZ, i)] += scalar * (valR * vzR - valL * vzL);
  Dest[Tidx(ENER, i)] += scalar * (valR - pR - (valL - pL));
}
