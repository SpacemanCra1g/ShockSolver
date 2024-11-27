#ifndef SR_CHARACTERISTICCLASS_H_
#define SR_CHARACTERISTICCLASS_H_
#include "../include/CharacteristicClass.hpp"

void Characteristics::EigenVectors(double *Prims, double *Cs, int Start,
                                   int Stop) {
  int i = 0;
  double vx, vy, vz, cs2, cs, v2tan, eta, d, h, lor, sroot;
  double **LL, **RR, *lam;

  for (i = Start; i < Stop; ++i) {
    vx = Prims[Tidx(VELX, i)];
    vy = Prims[Tidx(VELY, i)];
    vz = Prims[Tidx(VELZ, i)];
    d = Prims[Tidx(DENSP, i)];
    cs2 = Cs[i];
    LL = L[i];
    RR = R[i];
    lam = lambda[i];
    lor = 1.0 / std::sqrt(1.0 - (vx * vx + vy * vy + vz * vz));

#if EOS == IDEAL
    h = 1.0 + (GAMMA / (GAMMA - 1.0)) * Prims[Tidx(PRES, i)] / d;
#endif

    cs = std::sqrt(cs2);
    v2tan = vy * vy + vz * vz;
    eta = std::sqrt(1.0 - vx * vx - cs2 * v2tan);

    for (int ii = 0; ii < 5; ++ii) {
      for (int jj = 0; jj < 5; ++jj) {
        LL[ii][jj] = 0.0;
        RR[ii][jj] = 0.0;
      }
    }

    LL[0][1] = d * lor / (2.0 * cs * eta);
    LL[0][4] = 1.0 / (2.0 * h * cs2);

    LL[1][1] = -d * lor / (2.0 * cs * eta);
    LL[1][4] = 1.0 / (2.0 * h * cs2);

    LL[2][0] = 1.0;
    LL[2][4] = -1.0 / (h * cs2);

    LL[3][1] = vx * vy / (1.0 - vx * vx);
    LL[3][2] = 1.0;
    LL[3][4] = vy / (lor * lor * d * h * (1.0 - vx * vx));

    LL[4][1] = vx * vz / (1.0 - vx * vx);
    LL[4][3] = 1.0;
    LL[4][4] = vz / (lor * lor * d * h * (1.0 - vx * vx));

    RR[0][0] = 1.0;
    RR[0][1] = cs * eta / (d * lor);
    RR[0][2] =
        -cs * vy * (lor * eta * vx + cs) / (lor * lor * d * (1.0 - vx * vx));

    RR[0][3] =
        -cs * vz * (lor * eta * vx + cs) / (lor * lor * d * (1.0 - vx * vx));

    RR[0][4] = h * cs2;

    RR[1][0] = 1.0;
    RR[1][1] = -cs * eta / (d * lor);
    RR[1][2] =
        cs * vy * (lor * eta * vx - cs) / (lor * lor * d * (1.0 - vx * vx));

    RR[1][3] =
        cs * vz * (lor * eta * vx - cs) / (lor * lor * d * (1.0 - vx * vx));

    RR[1][4] = h * cs2;

    RR[2][0] = 1.0;
    RR[3][2] = 1.0;
    RR[4][3] = 1.0;

    sroot = cs2 / (lor * lor * (1 - cs2));
    lam[0] = (vx + std::sqrt(sroot * (1.0 - vx * vx + sroot))) / (1 + sroot);
    lam[1] = (vx - std::sqrt(sroot * (1.0 - vx * vx + sroot))) / (1 + sroot);
    lam[2] = vx;
    lam[3] = vx;
    lam[4] = vx;
  }
}

#endif // SR_CHARACTERISTICCLASS_H_
