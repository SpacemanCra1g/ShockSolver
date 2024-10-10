#include "../include/DomainClass.hpp"
#include "../include/SRVarConvert.hpp"

void Domain::Prims2Cons() {

  double P[NumVar], C[NumVar];

  for (int i = 0; i < xDim; ++i) {
    for (int var = 0; var < NumVar; ++var) {
      P[var] = Prims[Tidx(var, i)];
    }

    PrimConvert(P, C);

    for (int var = 0; var < NumVar; ++var) {
      Cons[Tidx(var, i)] = C[var];
    }
  }
}

void Domain::Cons2Prim() {

  double P[NumVar], C[NumVar];
  for (int i = 0; i < xDim; ++i) {
    for (int var = 0; var < NumVar; ++var) {
      C[var] = Cons[Tidx(var, i)];
    }

    ConConvert(C, P);

    for (int var = 0; var < NumVar; ++var) {
      Prims[Tidx(var, i)] = P[var];
    }
  }
}

void Domain::Press(int x) {
  double C[NumVar];
  for (int var = 0; var < NumVar; ++var) {
    C[var] = Cons[Tidx(var, x)];
  }
  PRES[x] = Pressure(C);
}

void Domain::SolvePressure() {
  for (int i = 0; i < xDim; ++i) {
    Press(i);
  }
}
