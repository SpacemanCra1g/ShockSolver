#include "../include/DomainClass.hpp"

void Domain::HLL_Speed(double *LP, double *RP, double *CsL, double *CsR, int i,
                       double &SL, double &SR) {
  double LLM, LRM, LLP, LRP;

  // LLM = LP[Tidx(VELX, i)] - std::sqrt(CsL[i]);
  // LLP = LP[Tidx(VELX, i)] + std::sqrt(CsL[i]);

  SignalSpeed(LP,CsL,i,LLM,LLP);
  SignalSpeed(RP,CsR,i+1,LRM,LRP);

  // LRM = RP[Tidx(VELX, i + 1)] - std::sqrt(CsR[i + 1]);
  // LRP = RP[Tidx(VELX, i + 1)] + std::sqrt(CsR[i + 1]);

  SL = std::fmin(LLM, LRM);
  SR = std::fmax(LLP, LRP);
}

void Domain::FillFlux(double *LP, double *FL, double *RP, double *FR, int i) {
  double d, vx, p, E;
  d = LP[Tidx(DENSP, i)];
  vx = LP[Tidx(VELX, i)];
  p = LP[Tidx(PRES, i)];

  FL[DENS] = d * vx;
  FL[MOMX] = d * vx * vx + p;
  E = d * (0.5 * vx * vx + p / ( (GAMMA - 1.0) * d));
  FL[ENER] = vx*(E + p);

  d = RP[Tidx(DENSP, i + 1)];
  vx = RP[Tidx(VELX, i + 1)];
  p = RP[Tidx(PRES, i + 1)];

  FR[DENS] = d * vx;
  FR[MOMX] = d * vx * vx + p;
  E = d * (0.5 * vx * vx + p / ( (GAMMA - 1.0) * d));
  FR[ENER] = vx*(E + p);
}
