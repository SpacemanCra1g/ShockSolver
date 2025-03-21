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
  double vx, cs2, ca, Bdot, Bx, By, Bz, d, p;

  vx = Uin[Tidx(VELX, i)];
  d = Uin[Tidx(DENS, i)];
  p = Uin[Tidx(PRES, i)];
  Bx = Uin[Tidx(BX, i)];
  By = Uin[Tidx(BY, i)];
  Bz = Uin[Tidx(BZ, i)];
  // cs2 = CS[i];

  Bdot = Bx * Bx + By * By + Bz * Bz;

  // ca = (cs2 + Bdot / d) + std::sqrt(std::pow(cs2 - Bdot / (d * d), 2) +
  //                                   4 * cs2 * (By * By + Bz * Bz) / d);
  // ca *= 0.5;
  // ca = std::sqrt(ca);
  cs2 = .5 * ((GAMMA * p + Bdot) / d +
              std::sqrt(std::pow((GAMMA * p + Bdot) / d, 2) -
                        (4 * GAMMA * p * Bx * Bx) / (d * d)));
  // cs2 = .5 *
  //       (GAMMA * p + Bdot +
  //        std::sqrt(std::pow(GAMMA * p - Bdot, 2) +
  //                  4 * GAMMA * p * (By * By + Bz * Bz))) /
  //       d;

  ca = std::sqrt(cs2);
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
