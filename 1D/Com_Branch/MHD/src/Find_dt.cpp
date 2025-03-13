#include "../include/DomainClass.hpp"
// #include "../include/SRVarConvert.hpp"

double FindMaximum(double *Array, const int size) {
  double Max = 1.e-10;
  for (int i = XStart - 1; i < size; ++i) {
    if (Max < Array[i]) {
      Max = Array[i];
    }
  }
  return Max;
}

void Domain::SignalSpeed(double *Uin, double *CS, int i, double &CSL,
                         double &CSR) {
  double vx, cs2, ca, cfP, cfM, Bmag;

  vx = Uin[Tidx(VELX, i)];
  cs2 = CS[i];

  Bmag = std::sqrt(Uin[Tidx(BX, i)] * Uin[Tidx(BX, i)] +
                   Uin[Tidx(BY, i)] * Uin[Tidx(BY, i)] +
                   Uin[Tidx(BZ, i)] * Uin[Tidx(BZ, i)]);

  ca = Uin[Tidx(BX, i)] / std::sqrt(Uin[Tidx(DENS, i)]);
  cfP = (GAMMA * Uin[Tidx(PRES, i)] + Bmag) / Uin[Tidx(DENS)];
  cfP = std::sqrt(
      .5 * (cfP + std::sqrt(cfP * cfP -
                            4 * GAMMA * Uin[Tidx(PRES, i)] * Uin[Tidx(BX, i)] *
                                Uin[Tidx(BX, i)] /
                                (Uin[Tidx(PRES, i)] * Uin[Tidx(PRES, i)]))));

  cfM = (GAMMA * Uin[Tidx(PRES, i)] + Bmag) / Uin[Tidx(DENS)];
  cfM = std::sqrt(
      .5 * (cfM + std::sqrt(cfM * cfM -
                            4 * GAMMA * Uin[Tidx(PRES, i)] * Uin[Tidx(BX, i)] *
                                Uin[Tidx(BX, i)] /
                                (Uin[Tidx(PRES, i)] * Uin[Tidx(PRES, i)]))));

  cfP = std::fmax(std::fabs(cfM), std::fabs(cfP));
  ca = std::fmax(std::fabs(cfP), std::fabs(ca));

  CSR = vx + ca;
  CSL = vx - ca;
}

void Domain::Find_Cs(double *Uin, double *CS, int start, int end) {

  for (int i = start; i < end; ++i) {
    CS[i] = GAMMA * Uin[Tidx(PRES, i)] / Uin[Tidx(DENS, i)];
  }
}

void Domain::Find_dt() {
  double CsL, CsR;
  Cons2Prim(Cons, Prims, XStart - 1, XEnd + 1);

  Find_Cs(Prims, Cs, XStart - 1, XEnd + 1);

  for (int i = XStart - 1; i < XEnd + 1; ++i) {

    SignalSpeed(Prims, Cs, i, CsL, CsR);

    Buffer[i] = std::fmax(std::fabs(CsL), std::fabs(CsR));
  }

  dt = dx / FindMaximum(Buffer, XEnd + 1);
  dt *= CFL;

  if (T + dt > TN) {
    dt = TN - T;
  }

#if SlowStart
  if (dt > dt_sim) {
    dt = dt_sim;
    dt_sim *= 2.0;
  }
#endif
}
