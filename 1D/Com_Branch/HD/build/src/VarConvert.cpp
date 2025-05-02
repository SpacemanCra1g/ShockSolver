#include "../include/DomainClass.hpp"
#include "../include/Parameters.h"


void Domain::Prims2Cons(double *Uin, double *Uout, int start, int stop) {
  double d, vx, p;
  for (int i = start; i < stop; ++i) {
    d = Uin[Tidx(DENSP, i)];
    vx = Uin[Tidx(VELX, i)];
    p = Uin[Tidx(PRES, i)];

    Uout[Tidx(DENS, i)] = d;
    Uout[Tidx(MOMX, i)] = d * vx;
    Uout[Tidx(ENER, i)] = d * (vx * vx * 0.5 + p / ((GAMMA - 1.0) * d));
  }
}

int Domain::Cons2Prim(double *Uin, double *Uout, int start, int stop) {
  double d, mx, E, vx;

  for (int i = start; i < stop; ++i) {
    d = Uin[Tidx(DENS, i)];
    mx = Uin[Tidx(MOMX, i)];
    E = Uin[Tidx(ENER, i)];

    Uout[Tidx(DENSP, i)] = d;
    vx = mx / d;
    Uout[Tidx(VELX, i)] = vx;
    Uout[Tidx(PRES, i)] = (E / d - vx * vx * 0.5) * (GAMMA - 1.0) * d;
    ConversionFailed[i] = false;
  }
  return 0;
}

