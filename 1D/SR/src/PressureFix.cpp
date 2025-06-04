#include "../include/DomainClass.hpp"

// This file is ripped almost completely from PLUTO/Src/RHD/rhd_pressure_fix.c
// So all credit goes to those authors

#define Min_p (1.e-2)
#define Max_Iter 20

int Domain::PressureFix(double *Uin, double *Uout, int i) {
  double D, m2, m, umax, u0, lor, plor, Dh, f0, f1, u1, alpha, du;
  const double p = Min_p;
  int done, k;

  D = Uin[Tidx(DENS, i)];
  m2 = Uin[Tidx(MOMX, i)] * Uin[Tidx(MOMX, i)] +
       Uin[Tidx(MOMY, i)] * Uin[Tidx(MOMY, i)] +
       Uin[Tidx(MOMZ, i)] * Uin[Tidx(MOMZ, i)];
  m = std::sqrt(m2);
  umax = m / D;
  u0 = umax;

  lor = std::sqrt(1.0 + u0 * u0);
  plor = p * lor;

#if EOS == IDEAL
  alpha = GAMMA / (GAMMA - 1.0);
  Dh = D + plor * alpha;
#endif

  f0 = m / Dh - u0;

  u1 = (-D + std::sqrt(D * D + 4.0 * m * alpha * p)) / (2.0 * alpha * p);

  done = 0;

  for (k = 1; k < Max_Iter; k++) {
    lor = sqrt(1.0 + u1 * u1);
    plor = p * lor;
#if EOS == IDEAL
    Dh = D + plor * alpha;
#elif EOS == TAUB
    Dh = 2.5 * plor + sqrt(2.25 * plor * plor + D * D);
#endif
    f1 = m / Dh - u1;

    if (done == 1)
      break;
    du = (u1 - u0) / (f1 - f0) * f1;
    u0 = u1;
    f0 = f1;

    u1 -= du;
    u1 = std::fmin(u1, umax);
    u1 = std::fmax(u1, 0.0);
    if (fabs(f1) < 1.e-9)
      done = 1;
  }

  if (k >= Max_Iter) {
    return 2;
  }

  lor = sqrt(1.0 + u1 * u1);
  plor = p * lor;
#if EOS == IDEAL
  Dh = D + plor * alpha;
#endif

  Uout[Tidx(DENSP, i)] = Uin[Tidx(DENS, i)] / lor;
  Uout[Tidx(PRES, i)] = p;
  Uin[Tidx(ENER, i)] =
      Dh * lor - p; // Redefining energy to match the fixed pressure

  f0 = 1.0 / (Dh * lor); /* = 1 / W */

  for (int var = MOMX; var <= MOMZ; ++var) {
    Uout[Tidx(var, i)] = Uin[Tidx(var, i)] * f0;
  }

  return 0;
}
#undef Max_Iter
#undef Min_p
