#include "../include/DomainClass.hpp"
#include "../include/SRVarConvert.hpp"

double FindMinimum(double *Array, const int size) {
  double Min = 1E100;
  for (int i = 0; i < size; ++i) {
    if (Min > Array[i]) {
      Min = Array[i];
    }
  }
  return Min;
}

void Domain::Find_Cs() {

  // for (int i = 0; i < xDim; ++i) {
  //   // Cs[i] = HR_CS(i);
  //   Cs[i] = SRHD_CS(i);
  //   std::cout << Cs[i] << std::endl;
  // }
  // exit(0);
}

void Domain::Find_dt() {
  Cons2Prim();
  // Find_Cs();
  double Cs;
  double CsL, CsR;
  double P[NumVar];

  for (int i = 0; i < xDim; ++i) {
    for (int var = 0; var < NumVar; ++var) {
      P[var] = Prims[Tidx(var, i)];
    }

    Cs = SRHD_CS(i);

    SignalSpeed(P, Cs, CsL, CsR);

    Buffer[i] = std::fmax(dx / CsL, dx / CsR);
  }

  dt = FindMinimum(Buffer, xDim);
  dt *= CFL;

  if (T + dt > TN) {
    dt = TN - T;
  }

#if SlowStart == True
  if (dt > dt_sim) {
    dt = dt_sim;
    dt_sim *= 2.0;
  }
#endif
}

double Domain::HD_CS(int i) { return std::sqrt(PRES[i] * GAMMA / DENS[i]); }

double Domain::SRHD_CS(int i) {
  double P[NumVar];
  for (int var = 0; var < NumVar; ++var) {
    P[var] = Prims[Tidx(var, i)];
  }
  return SRH_CS(P);
}
