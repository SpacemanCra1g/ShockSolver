#include "../include/DomainClass.hpp"
// #include "../include/SRVarConvert.hpp"

double FindMinimum(double *Array, const int size) {
  double Min = 1E100;
  for (int i = 0; i < size; ++i) {
    if (Min > Array[i]) {
      Min = Array[i];
    }
  }
  return Min;
}

void Domain::SignalSpeed(double *Uin, int i, double &CSL, double &CSR) {
  double v2, Delt2, Nu2, vx, vy, vz, sroot, cs2, lor;

  vx = Uin[Tidx(VELX, i)];
  vy = Uin[Tidx(VELY, i)];
  vz = Uin[Tidx(VELZ, i)];
  cs2 = Cs[i];

  v2 = vx * vx + vy * vy + vz * vz;

  lor = 1.0 / (std::sqrt(1.0 - v2));

  Delt2 = 1.0 - v2 * cs2;
  Nu2 = 1.0 - vx * vx - cs2 * (vy * vy + vz * vz);

  sroot = std::sqrt(cs2 * Nu2);

  CSR = (lor * vx * (1.0 - cs2) - sroot) / (lor * Delt2);
  CSL = (lor * vx * (1.0 - cs2) + sroot) / (lor * Delt2);
}

void Domain::Find_Cs(double *Uin, int start, int end) {
  double h;
  for (int i = start; i < end; ++i) {
#if EOS == IDEAL
    h = 1 + (GAMMA / (GAMMA - 1)) * Uin[Tidx(PRES, i)] / Uin[Tidx(DENSP, i)];
    Cs[i] = GAMMA * Uin[Tidx(PRES, i)] / (h * Uin[Tidx(DENSP, i)]);
#endif
  }
}

void Domain::Find_dt() {
  double CsL, CsR;
  Cons2Prim(Cons, Prims, 0, xDim);
  Find_Cs(Prims, 0, xDim);

  for (int i = 0; i < xDim; ++i) {

    SignalSpeed(Cs, i, CsL, CsR);

    Buffer[i] = std::fmax(dx / CsL, dx / CsR);
  }

  dt = FindMinimum(Buffer, xDim);
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

// double Domain::HD_CS(int i) { return std::sqrt(Pres[i] * GAMMA / Dens[i]); }

// double Domain::SRHD_CS(int i) {
//   double P[NumVar];
//   for (int var = 0; var < NumVar; ++var) {
//     P[var] = Prims[Tidx(var, i)];
//   }
//   return SRH_CS(P);
// }
