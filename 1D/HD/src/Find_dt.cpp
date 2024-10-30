#include "../include/DomainClass.hpp"

double FindMaximum(double *Array, const int size) {
  double Min = 1e-100;
  for (int i = XStart - 1; i < size; ++i) {

    Min = (Min < Array[i]) ? Array[i] : Min;
  }
  return Min;
}

void Domain::SignalSpeed(double *Uin, double *CS, int i, double &CSL,
                         double &CSR) {

  double u = Uin[Tidx(VELX, i)];
  CSR = u + CS[i];
  CSL = u - CS[i];
}

void Domain::Find_Cs(double *Uin, double *CS, int start, int end) {

  for (int i = start; i < end; ++i) {
#if EOS == IDEAL
    CS[i] = std::sqrt(Uin[Tidx(PRES, i)] * GAMMA / Uin[Tidx(DENSP, i)]);
#endif
  }
}

void Domain::Find_dt() {
  int err;
  err = Cons2Prim(Cons, Prims, 0, xDim);

  if (err) {
    std::cout << "Con Conversion failure in the find dt stage" << std::endl;
    exit(0);
  }
  Find_Cs(Prims, Cs, 0, xDim);

  for (int i = 0; i < xDim; ++i) {

    // SignalSpeed(Prims, Cs, i, CsL, CsR);

    // Buffer[i] = std::fmax(std::fabs(CsL), std::fabs(CsR));
    Buffer[i] = (std::fabs(Xvel[i]) + Cs[i]);
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
