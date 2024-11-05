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
  double vx, vy, vz, cs2, v2, sroot; // Delt2, Nu2,
  double lor;

  vx = Uin[Tidx(VELX, i)];
  vy = Uin[Tidx(VELY, i)];
  vz = Uin[Tidx(VELZ, i)];
  cs2 = CS[i];

  v2 = vx * vx + vy * vy + vz * vz;
  lor = 1.0 / std::sqrt(1.0 - v2);

  sroot = cs2 / (lor * lor * (1 - cs2));
  CSR = (vx + std::sqrt(sroot * (1.0 - vx * vx + sroot))) / (1 + sroot);
  CSL = (vx - std::sqrt(sroot * (1.0 - vx * vx + sroot))) / (1 + sroot);

  // The other attept
  // sroot =
  //     std::sqrt(cs2 * (1.0 - vx * vx - (vy + vy + vz * vz) * cs2 * (1.0 -
  //     v2)));
  // CSR = (vx * (1.0 - cs2) + sroot) / (1.0 - v2 * cs2);
  // CSL = (vx * (1.0 - cs2) - sroot) / (1.0 - v2 * cs2);
}

void Domain::Find_Cs(double *Uin, double *CS, int start, int end) {
  double h;
  for (int i = start; i < end; ++i) {
#if EOS == IDEAL
    h = 1.0 +
        (GAMMA / (GAMMA - 1.0)) * Uin[Tidx(PRES, i)] / Uin[Tidx(DENSP, i)];
    CS[i] = GAMMA * Uin[Tidx(PRES, i)] / (h * Uin[Tidx(DENSP, i)]);
#endif
  }
}

void Domain::Find_dt() {
  double CsL, CsR;
  int err;
  err = Cons2Prim(Cons, Prims, XStart - 1, XEnd + 1);
  if (err) {
    std::cout << "Con Conversion failure in the find dt stage" << std::endl;
    exit(0);
  }
  Find_Cs(Prims, Cs, XStart - 1, XEnd + 1);

  for (int i = XStart - 1; i < XEnd + 1; ++i) {

    SignalSpeed(Prims, Cs, i, CsL, CsR);

    Buffer[i] = std::fmax(std::fabs(CsL), std::fabs(CsR));
    Buffer[i] = std::fmax(Buffer[i], std::fabs(Prims[Tidx(VELX, i)]));
  }

  dt = dx / FindMaximum(Buffer, XEnd + 1);
  dt *= CFL * .8;

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
