#ifndef SR_CHARACTERISTICCLASS_H_
#define SR_CHARACTERISTICCLASS_H_
#include "../include/CharacteristicClass.hpp"

void Characteristics::EigenVectors(double *Prims, double *Cs, int Start,
                                   int Stop) {
  int i = 0;
  double vx, cs2, cs, d;
  double **LL, **RR, *lam;

  for (i = Start; i < Stop; ++i) {
    vx = Prims[Tidx(VELX, i)];
    d = Prims[Tidx(DENSP, i)];
    cs2 = Cs[i];
    LL = L[i];
    RR = R[i];
    lam = lambda[i];

    cs = std::sqrt(cs2);

    for (int ii = 0; ii < NumVar; ++ii) {
      for (int jj = 0; jj < NumVar; ++jj) {
        LL[ii][jj] = 0.0;
        RR[ii][jj] = 0.0;
      }
    }

    LL[0][1] = 1.0;
    LL[0][2] = -1.0 / (d * cs);

    LL[1][0] = 1.0;
    LL[1][2] = -1.0 / cs2;

    LL[2][1] = 1.0;
    LL[2][2] = 1.0 / (d * cs);

    RR[0][0] = -d / (2.0 * cs);
    RR[0][1] = 0.5;
    RR[0][2] = -d * cs * 0.5;

    RR[1][0] = 1.0;

    RR[2][0] = d / (2.0 * cs);
    RR[2][1] = 0.5;
    RR[2][2] = d * cs * 0.5;

    lam[0] = vx - cs;
    lam[1] = vx;
    lam[2] = vx + cs;
  }
}

#endif // SR_CHARACTERISTICCLASS_H_
