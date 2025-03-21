#include "../include/DomainClass.hpp"
#include "../include/Parameters.h"

void Domain::Prims2Cons(double *Uin, double *Uout, int start, int stop) {
  double d, vx, vy, vz, p, Bx, By, Bz, E, Z;
  for (int i = start; i < stop; ++i) {
    d = Uin[Tidx(DENSP, i)];
    vx = Uin[Tidx(VELX, i)];
    vy = Uin[Tidx(VELY, i)];
    vz = Uin[Tidx(VELZ, i)];
    p = Uin[Tidx(PRES, i)];
    Bx = Uin[Tidx(BX, i)];
    By = Uin[Tidx(BY, i)];
    Bz = Uin[Tidx(BZ, i)];

    Uout[Tidx(DENS, i)] = d;
    Uout[Tidx(MOMX, i)] = d * vx;
    Uout[Tidx(MOMY, i)] = d * vy;
    Uout[Tidx(MOMZ, i)] = d * vz;

    E = (vx * vx + vy * vy + vz * vz) * 0.5 + p / ((GAMMA - 1.0) * d);
    Z = E + .5 * (Bx * Bx + By * By + Bz * Bz) / d;
    Uout[Tidx(ENER, i)] = d * Z;

    Uout[Tidx(BX, i)] = Bx;
    Uout[Tidx(BY, i)] = By;
    Uout[Tidx(BZ, i)] = Bz;
  }
}

int Domain::Cons2Prim(double *Uin, double *Uout, int start, int stop) {
  double d, mx, my, mz, Bx, By, Bz, E, Z, vx, vy, vz;

  for (int i = start; i < stop; ++i) {
    d = Uin[Tidx(DENS, i)];
    mx = Uin[Tidx(MOMX, i)];
    my = Uin[Tidx(MOMY, i)];
    mz = Uin[Tidx(MOMZ, i)];
    Z = Uin[Tidx(ENER, i)];
    Bx = Uin[Tidx(BX, i)];
    By = Uin[Tidx(BY, i)];
    Bz = Uin[Tidx(BZ, i)];

    Uout[Tidx(BX, i)] = Bx;
    Uout[Tidx(BY, i)] = By;
    Uout[Tidx(BZ, i)] = Bz;
    Uout[Tidx(DENSP, i)] = d;

    vx = mx / d;
    vy = my / d;
    vz = mz / d;

    Uout[Tidx(VELX, i)] = vx;
    Uout[Tidx(VELY, i)] = vy;
    Uout[Tidx(VELZ, i)] = vz;
    E = (Z - .5 * (Bx * Bx + By * By + Bz * Bz)) / d;
    Uout[Tidx(PRES, i)] =
        (E - (vx * vx + vy * vy + vz * vz) * 0.5) * (GAMMA - 1.0) * d;

    ConversionFailed[i] = false;
  }
  return 0;
}

// void Domain::Press(int x) {
//   double C[NumVar];
//   for (int var = 0; var < NumVar; ++var) {
//     C[var] = Cons[Tidx(var, x)];
//   }
//   // PRES[x] = Pressure(C);n
// }

// void Domain::SolvePressure() {
//   for (int i = 0; i < xDim; n++ i) {
//     Press(i);
//   }
// }
