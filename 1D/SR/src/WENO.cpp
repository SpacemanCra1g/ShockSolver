#include "../include/DomainClass.hpp"

void Domain::Weno(int xdir) {

  // Probably slow to do this pointer aliasing, but certainly makes
  // the code much easier to write
  for (int var = 0; var < NumVar; ++var) {

    double *Center = &Cons[Tidx(var, xdir)];

    double p1L = (-1.0 / 6.0) * (Center[-2]) + (5.0 / 6.0) * (Center[-1]) +
                 (1.0 / 3.0) * (Center[0]);

    double p1R = (1.0 / 3.0) * (Center[-2]) + (-7.0 / 6.0) * (Center[-1]) +
                 (11.0 / 6.0) * (Center[0]);

    double p2L = (1.0 / 3.0) * (Center[-1]) + (5.0 / 6.0) * (Center[0]) +
                 (-1.0 / 6.0) * (Center[1]);

    double p2R = (-1.0 / 6.0) * (Center[-1]) + (5.0 / 6.0) * (Center[0]) +
                 (1.0 / 3.0) * (Center[1]);

    double p3L = (11.0 / 6.0) * (Center[0]) + (-7.0 / 6.0) * (Center[1]) +
                 (1.0 / 3.0) * (Center[2]);

    double p3R = (1.0 / 3.0) * (Center[0]) + (5.0 / 6.0) * (Center[1]) +
                 (-1.0 / 6.0) * (Center[2]);

    double Beta1 =
        (13.0 / 12.0) *
            (std::pow(Center[-2] - 2.0 * Center[-1] + Center[0], 2)) +
        0.25 * (std::pow(Center[-2] - 4.0 * Center[-1] + 3.0 * Center[0], 2));

    double Beta2 = (13.0 / 12.0) *
                       (std::pow(Center[-1] - 2.0 * Center[0] + Center[1], 2)) +
                   0.25 * (std::pow(Center[-1] - Center[1], 2));

    double Beta3 =
        (13.0 / 12.0) * (std::pow(Center[0] - 2.0 * Center[1] + Center[2], 2)) +
        0.25 * (std::pow(3.0 * Center[0] - 4.0 * Center[1] + Center[2], 2));

    double eps = 1E-36;

    double w1L, w2L, w3L, w1R, w2R, w3R, wLSum, wRSum;

    w1L = 0.3 / (eps + Beta1);
    w1R = 0.1 / (eps + Beta1);

    w2L = 0.6 / (eps + Beta2);
    w2R = 0.6 / (eps + Beta2);

    w3L = 0.1 / (eps + Beta3);
    w3R = 0.3 / (eps + Beta3);

    wLSum = w1L + w2L + w3L;
    wRSum = w1R + w2R + w3R;

    w1L /= wLSum;
    w1R /= wRSum;

    w2L /= wLSum;
    w2R /= wRSum;

    w3L /= wLSum;
    w3R /= wRSum;

    // Populate the Fluxes, of each variable

    FluxWalls_Cons[LEFT][Tidx(var, xdir)] = w1L * p1L + w2L * p2L + w3L * p3L;
    FluxWalls_Cons[RIGHT][Tidx(var, xdir)] = w1R * p1R + w2R * p2R + w3R * p3R;

    // if (var == PRES) {
    //   double PresTest =
    //       (Cons[Tidx(Pres, xdir)] - FluxDir[Left][qp][Pres][idx(xdir)]) *
    //       (FluxDir[Right][qp][Pres][idx(xdir)] - Cons[Tidx(Pres, xdir)]);

    //   double DensTest =
    //       (Cons[Tidx(Dens, xdir)] - FluxDir[Left][qp][Dens][idx(xdir)]) *
    //       (FluxDir[Right][qp][Dens][idx(xdir)] - Cons[Tidx(Dens, xdir)]);

    //   if (PresTepst < 0.0 || DensTest < 0.0) {
    //     for (int var = Dens; var <= Ener; ++var) {
    //       FluxDir[Left][qp][var][idx(xdir)] = Cons[Tidx(var, xdir)];
    //       FluxDir[Right][qp][var][idx(xdir)] = Cons[Tidx(var, xdir)];
    //     }
    //   }
    // }
  }
}
