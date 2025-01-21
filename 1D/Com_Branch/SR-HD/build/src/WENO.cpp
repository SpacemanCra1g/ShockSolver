#include "../include/DomainClass.hpp"

void Domain::Weno(int start, int stop) {
  double *Center, p1L, p1R, p2L, p2R, p3L, p3R, Beta1, Beta2, Beta3, eps, w1L,
      w2L, w3L, w1R, w2R, w3R, wLSum, wRSum, PresTest, DensTest;
  double **Stencil;
  double LeftState[5], RightState[5];

  Stencil = new double *[NumVar];
  for (int i = 0; i < 5; ++i) {
    Stencil[i] = new double[5];
  }

  // Probably slow to do this pointer aliasing, but certainly makes
  // the code much easier to write

#if EvolveChars
  Find_Cs(Prims, Cs, 0, REdgeX);
  Chars.EigenVectors(Prims, Cs, 0, REdgeX);
#endif

  for (int xdir = start; xdir < stop; ++xdir) {
    for (int var = 0; var < NumVar; ++var) {
      for (int stencil = 0; stencil < 5; ++stencil) {
        Stencil[var][stencil] = Prims[Tidx(var, xdir + stencil - 2)];
      }
    }

#if EvolveChars
    Chars.Char_Stencil(Stencil, xdir, 2);
#endif

    for (int var = 0; var < NumVar; ++var) {

      Center = &Stencil[var][2];

      p1L = (-1.0 / 6.0) * (Center[-2]) + (5.0 / 6.0) * (Center[-1]) +
            (1.0 / 3.0) * (Center[0]);

      p1R = (1.0 / 3.0) * (Center[-2]) + (-7.0 / 6.0) * (Center[-1]) +
            (11.0 / 6.0) * (Center[0]);

      p2L = (1.0 / 3.0) * (Center[-1]) + (5.0 / 6.0) * (Center[0]) +
            (-1.0 / 6.0) * (Center[1]);

      p2R = (-1.0 / 6.0) * (Center[-1]) + (5.0 / 6.0) * (Center[0]) +
            (1.0 / 3.0) * (Center[1]);

      p3L = (11.0 / 6.0) * (Center[0]) + (-7.0 / 6.0) * (Center[1]) +
            (1.0 / 3.0) * (Center[2]);

      p3R = (1.0 / 3.0) * (Center[0]) + (5.0 / 6.0) * (Center[1]) +
            (-1.0 / 6.0) * (Center[2]);

      Beta1 =
          (13.0 / 12.0) *
              (std::pow(Center[-2] - 2.0 * Center[-1] + Center[0], 2)) +
          0.25 * (std::pow(Center[-2] - 4.0 * Center[-1] + 3.0 * Center[0], 2));

      Beta2 = (13.0 / 12.0) *
                  (std::pow(Center[-1] - 2.0 * Center[0] + Center[1], 2)) +
              0.25 * (std::pow(Center[-1] - Center[1], 2));

      Beta3 =
          (13.0 / 12.0) *
              (std::pow(Center[0] - 2.0 * Center[1] + Center[2], 2)) +
          0.25 * (std::pow(3.0 * Center[0] - 4.0 * Center[1] + Center[2], 2));

      eps = 1E-36;

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

      LeftState[var] = w1L * p1L + w2L * p2L + w3L * p3L;
      RightState[var] = w1R * p1R + w2R * p2R + w3R * p3R;
    }
#if EvolveChars

    Chars.Char_Recover_Prims(LeftState, FluxWalls_Prims[LEFT], xdir);
    Chars.Char_Recover_Prims(RightState, FluxWalls_Prims[RIGHT], xdir);

#endif
    PresTest = (Prims[Tidx(PRES, xdir)] - LeftState[PRES]) *
               (RightState[PRES] - Prims[Tidx(PRES, xdir)]);

    DensTest = (Prims[Tidx(DENS, xdir)] - LeftState[DENS]) *
               (RightState[DENS] - Prims[Tidx(DENS, xdir)]);

    if (PresTest < 0.0 || DensTest < 0.0) {
      for (int var = 0; var < NumVar; ++var) {
        FluxWalls_Prims[LEFT][Tidx(var, xdir)] = Prims[Tidx(var, xdir)];
        FluxWalls_Prims[RIGHT][Tidx(var, xdir)] = Prims[Tidx(var, xdir)];
      }
    } else {

      for (int var = 0; var < NumVar; ++var) {
        FluxWalls_Prims[LEFT][Tidx(var, xdir)] = LeftState[var];
        FluxWalls_Prims[RIGHT][Tidx(var, xdir)] = RightState[var];
      }
    }
  }
}
