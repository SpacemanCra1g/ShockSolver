#include "../include/DomainClass.hpp"
#include "../include/Parameters.h"

void Domain::Prims2Cons(double *Uin, double *Uout, int start, int stop) {
  double vx, p, d, e;
  for (int i = start; i < stop; ++i) {
    d = Uin[Tidx(DENSP, i)];
    p = Uin[Tidx(PRES, i)];
    vx = Uin[Tidx(VELX, i)];

#if EOS == IDEAL
    e = p / ((GAMMA - 1.0) * d);
#endif

    Uout[Tidx(DENS, i)] = d;
    Uout[Tidx(MOMX, i)] = d * vx;
    Uout[Tidx(ENER, i)] = (vx * vx * d * 0.5) + e * d;
  }
}

int Domain::Cons2Prim(double *Uin, double *Uout, int start, int stop) {
  double mx, d, E, vx, e;
  for (int i = start; i < stop; ++i) {
    E = Uin[Tidx(ENER, i)];
    d = Uin[Tidx(DENS, i)];
    mx = Uin[Tidx(MOMX, i)];

    Uout[Tidx(DENSP, i)] = d;
    vx = mx / d;
    Uout[Tidx(VELX, i)] = vx;
#if EOS == IDEAL
    e = (E / d) - 0.5 * vx * vx;
#endif
    Uout[Tidx(PRES, i)] = (GAMMA - 1.0) * d * e;
  }
  return 0;
}
