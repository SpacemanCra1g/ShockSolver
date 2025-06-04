#include "../include/DomainClass.hpp"

#define MAX_ITER 20
#define TOL (1.e-9)

int Domain::EnergyInverter(double *Uin, double *Uout, int i) {

  double D, E, m2, Q, gam, m, p, D1, eps2, pmin;
  double alpha, alpha2, h, lor2, lor, tau, theta;
  double dh_dp, dh_dtau, dp, yp, dyp;
  int Idx;

  D = Uin[Tidx(DENS, i)];
  E = Uin[Tidx(ENER, i)];
  m2 = 0.0;
  for (int var = MOMX; var <= MOMZ; ++var) {
    m2 += Uin[Tidx(var, i)] * Uin[Tidx(var, i)];
  }

  Q = E - std::sqrt(m2 + D * D);

  if (Q < 0.0) {
    // The equation does not admit a solution
    return 1;
  }

#if EOS == IDEAL
  gam = GAMMA / (GAMMA - 1);
#endif
  m = std::sqrt(m2);
  p = m - E;
  p = std::fmax(p, 0.0);
  D1 = 1.0 / D;

  eps2 = 1.e-12;
  pmin = std::sqrt(m2 / (1.0 - eps2)) - E;

  // We solve f(p) = 0 by Newton's method

  p = std::fmax(p, pmin);

  for (Idx = 0; Idx < MAX_ITER; ++Idx) {
    alpha = E + p;
    alpha2 = alpha * alpha;
    lor2 = 1.0 - m2 / alpha2;
    lor2 = 1.0 / lor2;

    if (lor2 < 1.0) {
      std::cout << "ERROR In the LOR Factor" << std::endl;
      std::cout << "Lor2 = " << lor2 << std::endl;
      exit(0);
    }

    lor = std::sqrt(lor2);
    tau = lor * D1;
    theta = p * tau;

#if EOS == IDEAL
    h = 1.0 + gam * theta;
    dh_dp = gam * tau;
    dh_dtau = gam * p;
#endif
    yp = D * h * lor - E - p;
    dyp = D * lor * dh_dp -
          m2 * lor2 * lor / (alpha2 * alpha) * (lor * dh_dtau + D * h) - 1.0;
    dp = yp / dyp;
    p -= dp;

    if (p < pmin) {
      p = pmin;
    }
    if (std::fabs(dp) < TOL * E) {
      break;
    }
  }

  // Check Solution
  if (p < 0.0) {
    return 2;
  }
  if (std::isnan(p)) {
    return 4;
  }
  if (Idx >= MAX_ITER || std::fabs(yp / (E + p)) > 1.e-4 || p < (m - E)) {
    return 3;
  }

  // If accurate solution has been found, we update the vars
  alpha = E + p;
  alpha2 = alpha * alpha;
  lor2 = alpha2 / (alpha2 - m2);
  lor = std::sqrt(lor2);
  tau = lor * D1;
  theta = p * tau;

#if EOS == IDEAL
  h = 1.0 + gam * theta;
#endif

  Uout[Tidx(DENSP, i)] = 1.0 / tau;
  Uout[Tidx(PRES, i)] = p;
  double temp = 1.0 / (E + p);

  Uout[Tidx(VELX, i)] = Uin[Tidx(MOMX, i)] * temp;
  Uout[Tidx(VELY, i)] = Uin[Tidx(MOMY, i)] * temp;
  Uout[Tidx(VELZ, i)] = Uin[Tidx(MOMZ, i)] * temp;

  return 0;
}

#undef MAX_ITER
#undef TOL
