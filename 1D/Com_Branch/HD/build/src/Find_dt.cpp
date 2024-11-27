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
  double vx, cs2;

  vx = Uin[Tidx(VELX, i)];
  cs2 = CS[i];

  CSR = vx + std::sqrt(cs2);
  CSL = vx - std::sqrt(cs2);
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
