#include "../include/DomainClass.hpp"

void Domain::HLL_Speed(double *LP, double *RP, double *CsL, double *CsR, int i,
                       double &SL, double &SR) {

  // We will attempt the Davis Estimate
  double LLM, LRM, LLP, LRP;
  SignalSpeed(LP, CsL, i, LLM, LLP);
  SignalSpeed(RP, CsR, i + 1, LRM, LRP);

  SL = std::fmin(LLM, LRM);
  SR = std::fmax(LLP, LRP);
}

void Domain::FillFlux(double *LP, double *FL, double *RP, double *FR, int i) {
  double vx, vy, vz, d, p, Bx, By, Bz, TotP, E, Z;
  d = LP[Tidx(DENSP, i)];
  vx = LP[Tidx(VELX, i)];
  vy = LP[Tidx(VELY, i)];
  vz = LP[Tidx(VELZ, i)];
  p = LP[Tidx(PRES, i)];
  Bx = LP[Tidx(BX, i)];
  By = LP[Tidx(BY, i)];
  Bz = LP[Tidx(BZ, i)];
  TotP = p + .5 * (Bx * Bx + By * By + Bz * Bz) / d;

  E = p / ((GAMMA - 1) * d) + .5 * (vx * vx + vy * vy + vz * vz);
  Z = E + .5 * (Bx * Bx + By * By + Bz * Bz) / d;

  FL[DENS] = d * vx;
  FL[MOMX] = d * vx * vx + TotP - Bx * Bx;
  FL[MOMY] = d * vx * vy - Bx * By;
  FL[MOMZ] = d * vx * vz - Bx * Bz;
  FL[ENER] = vx * (d * Z + TotP) - Bx * (Bx * vx + vy * By + vz * Bz);
  FL[BX] = 0;
  FL[BY] = vx * By - Bx * vy;
  FL[BZ] = vx * Bz - Bx * vz;

  d = RP[Tidx(DENSP, i + 1)];
  vx = RP[Tidx(VELX, i + 1)];
  vy = RP[Tidx(VELY, i + 1)];
  vz = RP[Tidx(VELZ, i + 1)];
  p = RP[Tidx(PRES, i + 1)];
  Bx = RP[Tidx(BX, i + 1)];
  By = RP[Tidx(BY, i + 1)];
  Bz = RP[Tidx(BZ, i + 1)];
  TotP = p + .5 * (Bx * Bx + By * By + Bz * Bz) / d;

  E = p / ((GAMMA - 1) * d) + .5 * (vx * vx + vy * vy + vz * vz);
  Z = E + .5 * (Bx * Bx + By * By + Bz * Bz) / d;

  FR[DENS] = d * vx;
  FR[MOMX] = d * vx * vx + TotP - Bx * Bx;
  FR[MOMY] = d * vx * vy - Bx * By;
  FR[MOMZ] = d * vx * vz - Bx * Bz;
  FR[ENER] = vx * (d * Z + TotP) - Bx * (Bx * vx + vy * By + vz * Bz);
  FR[BX] = 0;
  FR[BY] = vx * By - Bx * vy;
  FR[BZ] = vx * Bz - Bx * vz;
}
