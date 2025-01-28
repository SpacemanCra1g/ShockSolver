#include "../include/DomainClass.hpp"
void Domain::Ausm(int Start, int Stop) {

  double rL, uL, pL, aL, ML, HL;
  double rR, uR, pR, aR, MR, HR;
  double Mp, Pp;
  double Mm, Pm;
  double Fp[3], Fm[3];
  double PosF, NegF;

  for (int i = Start; i < Stop; ++i) {
    rL = FluxWalls_Prims[RIGHT][Tidx(DENSP, i)];
    uL = FluxWalls_Prims[RIGHT][Tidx(VELX, i)];
    pL = FluxWalls_Prims[RIGHT][Tidx(PRES, i)];
    pL = (pL < 0.0) ? 0.001 : pL;
    aL = std::sqrt(GAMMA * pL / rL);
    ML = uL / aL;
    HL = aL * aL / (GAMMA - 1.0) + 0.5 * uL * uL;

    rR = FluxWalls_Prims[LEFT][Tidx(DENSP, i + 1)];
    uR = FluxWalls_Prims[LEFT][Tidx(VELX, i + 1)];
    pR = FluxWalls_Prims[LEFT][Tidx(PRES, i + 1)];
    pR = (pR < 0.0) ? 0.001 : pR;
    aR = std::sqrt(GAMMA * pR / rR);
    MR = uR / aR;
    HR = aR * aR / (GAMMA - 1.0) + 0.5 * uR * uR;

    if (ML <= -1.0) {
      Mp = 0.0;
      Pp = 0.0;
    } else if (ML < 1.0) {
      Mp = (ML + 1.0) * (ML + 1.0) * 0.25;
      Pp = pL * (1.0 + ML) * (1.0 + ML) * (2.0 - ML) * 0.25;
    } else {
      Mp = ML;
      Pp = pL;
    }

    if (MR <= -1.0) {
      Mm = MR;
      Pm = pR;
    } else if (MR < 1.0) {
      Mm = -(MR - 1.0) * (MR - 1.0) * 0.25;
      Pm = pR * (1.0 - MR) * (1.0 - MR) * (2.0 + MR) * 0.25;
    } else {
      Mm = 0.0;
      Pm = 0.0;
    }

    PosF = std::fmax(0.0, Mp + Mm) * aL;
    NegF = std::fmin(0.0, Mp + Mm) * aR;

    Fp[0] = PosF * rL;
    Fp[1] = PosF * rL * uL + Pp;
    Fp[2] = PosF * rL * HL;

    Fm[0] = NegF * rR;
    Fm[1] = NegF * rR * uR + Pm;
    Fm[2] = NegF * rR * HR;

    for (int var = 0; var < NumVar; ++var) {
      CellFlux[Tidx(var, i)] = Fp[var] + Fm[var];
    }
  }
};

double sign(double M) {
  if (std::fabs(M) < 1.e-13) {
    return 0;
  } else {
    return M / std::fabs(M);
  }
}

// void Domain::Autsm(int Start, int Stop) {

//   double rL, uL, pL, aL, ML, HL, aStrL;
//   double rR, uR, pR, aR, MR, HR, aStrR;
//   double Mp, Pp;
//   double Mm, Pm;
//   double Fp[3], Fm[3];
//   double PosF, NegF;
//   double A;
//   double alpha = 3.0 / 16.0;
//   double beta = 1.0 / 8.0;

//   for (int i = Start; i < Stop; ++i) {
//     rL = FluxWalls_Prims[RIGHT][Tidx(DENSP, i)];
//     uL = FluxWalls_Prims[RIGHT][Tidx(VELX, i)];
//     pL = FluxWalls_Prims[RIGHT][Tidx(PRES, i)];
//     pL = (pL < 0.0) ? 0.001 : pL;
//     aL = std::sqrt(GAMMA * pL / rL);
//     HL = aL * aL / (GAMMA - 1.0) + 0.5 * uL * uL;
//     aStrL = std::sqrt(2.0 * (GAMMA - 1.0) * HL / (GAMMA + 1.0));

//     rR = FluxWalls_Prims[LEFT][Tidx(DENSP, i + 1)];
//     uR = FluxWalls_Prims[LEFT][Tidx(VELX, i + 1)];
//     pR = FluxWalls_Prims[LEFT][Tidx(PRES, i + 1)];
//     pR = (pR < 0.0) ? 0.001 : pR;
//     aR = std::sqrt(GAMMA * pR / rR);
//     HR = aR * aR / (GAMMA - 1.0) + 0.5 * uR * uR;
//     aStrR = std::sqrt(2.0 * (GAMMA - 1.0) * HR / (GAMMA + 1.0));

//     A = std::fmin(aStrL * aStrL / std::fmax(aStrL, std::fabs(uL)),
//                   aStrR * aStrR / std::fmax(aStrR, std::fabs(uR)));

//     ML = uL / A;
//     MR = uR / A;

//     if (std::fabs(ML) >= 1.0) {
//       Mp = 0.5 * (ML + std::fabs(ML));
//       Pp = .5 * (1.0 + sign(ML));
//     } else {
//       Mp = 0.25 * (ML + 1.0) * (ML + 1.0) +
//            beta * (ML * ML - 1.0) * (ML * ML - 1.0); // Paper calls this 1/2
//       Pp = .25 * (ML + 1.0) * (ML + 1.0) * (2.0 - ML) +
//            alpha * ML * (ML * ML - 1.0) * (ML * ML - 1.0);
//     }

//     if (std::fabs(MR) >= 1.0) {
//       Mm = .5 * (MR - std::fabs(MR));
//       Pm = .5 * (1.0 - sign(MR));
//     } else {
//       Mm = -0.25 * (MR - 1.0) * (MR - 1.0) -
//            beta * (MR * MR - 1.0) * (MR * MR - 1.0); // Paper calls this 1/2
//       Pm = .25 * (MR - 1.0) * (MR - 1.0) * (2.0 + MR) -
//            alpha * MR * (MR * MR - 1.0) * (MR * MR - 1.0);
//     }

//     PosF = A * .5 * (Mp + Mm + std::fabs(Mp + Mm));
//     NegF = A * .5 * (Mp + Mm - std::fabs(Mp + Mm));

//     Fp[0] = PosF * rL;
//     Fp[1] = PosF * rL * uL + Pp * pL;
//     Fp[2] = PosF * HL * rL;

//     Fm[0] = NegF * rR;
//     Fm[1] = NegF * rR * uR + Pm * pR;
//     Fm[2] = NegF * HR * rR;

//     for (int var = 0; var < NumVar; ++var) {
//       CellFlux[Tidx(var, i)] = Fp[var] + Fm[var];
//     }
//   }
// };

void Domain::Autsm(int Start, int Stop) {

  int i;

  double aL, ML, MpL, PpL, asL2, asL, atL;
  double aR, MR, MmR, PmR, asR2, asR, atR;
  double a, m, mp, mm;
  double rhoL, pL, uL, rhoR, pR, uR;
  double HL, HR, Fp[3], Fm[3];
  double alpha = 3.0 / 16.0, beta = 0.125;
  // double alpha = 0.0, beta = 0.0;

  for (i = Start; i < Stop; i++) {

    rhoL = FluxWalls_Prims[RIGHT][Tidx(DENS, i)];
    uL = FluxWalls_Prims[RIGHT][Tidx(VELX, i)];
    pL = FluxWalls_Prims[RIGHT][Tidx(PRES, i)];
    pL = (pL < 0.0) ? 0.001 : pL;

    rhoR = FluxWalls_Prims[LEFT][Tidx(DENS, i + 1)];
    uR = FluxWalls_Prims[LEFT][Tidx(VELX, i + 1)];
    pR = FluxWalls_Prims[LEFT][Tidx(PRES, i + 1)];
    pR = (pR < 0.0) ? 0.001 : pR;

    aL = std::sqrt(GAMMA * pL / rhoL);
    aR = std::sqrt(GAMMA * pR / rhoR);

    asL2 = uL * uL;
    asL2 = aL * aL / (GAMMA - 1.0) + 0.5 * asL2;
    asL2 *= 2.0 * (GAMMA - 1.0) / (GAMMA + 1.0);

    asR2 = uR * uR;
    asR2 = aR * aR / (GAMMA - 1.0) + 0.5 * asR2;
    asR2 *= 2.0 * (GAMMA - 1.0) / (GAMMA + 1.0);

    asL = std::sqrt(asL2);
    asR = std::sqrt(asR2);

    atL = asL2 / std::fmax(asL, std::fabs(uL));
    atR = asR2 / std::fmax(asR, std::fabs(uR));

    a = std::fmin(atL, atR);
    /*
        a = 0.5*(aL + aR);
    */
    /* --------------------------------------------
            define split Mach numbers
            define pressure terms
       -------------------------------------------- */

    ML = uL / a;
    if (std::fabs(ML) >= 1.0) {
      MpL = 0.5 * (ML + std::fabs(ML));
      PpL = ML > 0.0 ? 1.0 : 0.0;
    } else {
      MpL = 0.25 * (ML + 1.0) * (ML + 1.0) +
            beta * (ML * ML - 1.0) * (ML * ML - 1.0);
      PpL = 0.25 * (ML + 1.0) * (ML + 1.0) * (2.0 - ML) +
            alpha * ML * (ML * ML - 1.0) * (ML * ML - 1.0);
    }

    MR = uR / a;
    if (std::fabs(MR) >= 1.0) {
      MmR = 0.5 * (MR - std::fabs(MR));
      PmR = MR > 0.0 ? 0.0 : 1.0;
    } else {
      MmR = -0.25 * (MR - 1.0) * (MR - 1.0) -
            beta * (MR * MR - 1.0) * (MR * MR - 1.0);
      PmR = 0.25 * (MR - 1.0) * (MR - 1.0) * (2.0 + MR) -
            alpha * MR * (MR * MR - 1.0) * (MR * MR - 1.0);
    }

    m = MpL + MmR;

    mp = a * 0.5 * (m + std::fabs(m));
    mm = a * 0.5 * (m - std::fabs(m));

    /* -------------------------------------------------------------
                       Compute fluxes
       ------------------------------------------------------------- */

    // PosF = A * .5 * (mp + mm + std::fabs(mp + mm));
    // NegF = A * .5 * (mp + mm - std::fabs(mp + mm));

    // Redefining this as energy
    HR = rhoR * (0.5 * uR * uR + pR / ((GAMMA - 1.0) * rhoR));
    HL = rhoL * (0.5 * uL * uL + pL / ((GAMMA - 1.0) * rhoL));

    Fp[0] = mp * rhoL;
    Fp[1] = mp * rhoL * uL;
    Fp[2] = mp * (HL + pL);

    Fm[0] = mm * rhoR;
    Fm[1] = mm * rhoR * uR;
    Fm[2] = mm * (HR + pR);

    for (int var = 0; var < NumVar; ++var) {
      CellFlux[Tidx(var, i)] = Fp[var] + Fm[var];
    }
    CellFlux[Tidx(MOMX, i)] += PpL * pL + PmR * pR;
  }
};

void Domain::Autsmup(int Start, int Stop) {
  // This is now the Ausm+up solver

  double alpha = 3.0 / 16.0, beta = 0.125, Ku = 0.75, Kp = 0.25, sigma = 1.0;
  double rhoL, uL, pL, momL, EnL;
  double rhoR, uR, pR, momR, EnR;
  double PhiL[3], PhiR[3];
  double aL, aR, asL2, asR2, asL, asR, atL, atR, a;
  double Mm2, Mo2, Mo, fa;
  double mm2, mo2, ML, MR, Mp1L, Mm1R, Mp2L, Mm2L, Mp2R, Mm2R, Mp4L, Pp5L, Mm4R;
  double Pm5R, scrh, M, m;

  // First we do the Prims to Cons conversion of the Flux walls
  Prims2Cons(FluxWalls_Prims[LEFT], FluxWalls_Cons[LEFT], Start, Stop + 1);
  Prims2Cons(FluxWalls_Prims[RIGHT], FluxWalls_Cons[RIGHT], Start, Stop);

  for (int i = Start; i < Stop; ++i) {

    // Load Data for the Riemann problem
    rhoL = FluxWalls_Prims[LEFT][Tidx(DENS, i)];
    uL = FluxWalls_Prims[LEFT][Tidx(VELX, i)];
    pL = FluxWalls_Prims[LEFT][Tidx(PRES, i)];
    pL = (pL > 0.0) ? pL : 0.1;

    rhoR = FluxWalls_Prims[RIGHT][Tidx(DENS, i + 1)];
    uR = FluxWalls_Prims[RIGHT][Tidx(VELX, i + 1)];
    pR = FluxWalls_Prims[RIGHT][Tidx(PRES, i + 1)];
    pR = (pR > 0.0) ? pR : 0.1;

    momL = FluxWalls_Cons[LEFT][Tidx(MOMX, i)];
    EnL = FluxWalls_Cons[LEFT][Tidx(ENER, i)];

    momR = FluxWalls_Cons[RIGHT][Tidx(MOMX, i + 1)];
    EnR = FluxWalls_Cons[RIGHT][Tidx(ENER, i + 1)];

    // Convective flux vectors Phi
    PhiL[DENS] = 1.0;
    PhiL[MOMX] = uL;
    PhiL[ENER] = (EnL - pL) / rhoL;

    PhiR[DENS] = 1.0;
    PhiR[MOMX] = uR;
    PhiR[ENER] = (EnR - pR) / rhoR;

    // Calculate numerical sound speed
    aL = std::sqrt(GAMMA * pL / rhoL);
    aR = std::sqrt(GAMMA * pR / rhoR);

    asL2 = aL * aL / (GAMMA - 1.0) + 0.5 * uL * uL;
    asL2 *= 2.0 * (GAMMA - 1.0) / (GAMMA + 1.0);

    asR2 = aR * aR / (GAMMA - 1.0) + 0.5 * uR * uR;
    asR2 *= 2.0 * (GAMMA - 1.0) / (GAMMA + 1.0);

    asL = std::sqrt(asL2);
    asR = std::sqrt(asR2);

    atL = asL2 / std::fmax(asL, std::fabs(uL));
    atR = asR2 / std::fmax(asR, std::fabs(uR));

    a = std::fmin(atL, atR);

    // alpha & beta coefficients
    Mm2 = .5 * (uL * uL + uR * uR) / a /
          a; // This is highly questionable, it might be valid but its odd
    Mo2 = std::fmax(mm2, .4);
    Mo2 = std::fmin(1.0, Mo2);
    Mo = std::sqrt(mo2);
    fa = Mo * (2.0 - Mo);

    alpha = (3.0 / 16.0) * (-4.0 + 5.0 * fa * fa);
    beta = 1.0 / 8.0;

    // Split Mach Number, Split Pressure
    ML = uL / a;
    MR = uR / a;

    Mp1L = .5 * (ML + std::fabs(ML));
    Mm1R = .5 * (MR - std::fabs(MR));

    Mp2L = .25 * (ML + 1.0) * (ML + 1.0);
    Mm2L = -.25 * (ML - 1.0) * (ML - 1.0);

    Mp2R = .25 * (MR + 1.0) * (MR + 1.0);
    Mm2R = -.25 * (MR - 1.0) * (MR - 1.0);

    if (std::fabs(ML) >= 1.0) {
      Mp4L = Mp1L;
      Pp5L = Mp1L / ML;
    } else {
      Mp4L = Mp2L * (1.0 - 16.0 * beta * Mm2L);
      Pp5L = Mp2L * ((2.0 - ML) - 16.0 * alpha * ML * Mm2L);
    }

    if (std::fabs(MR) >= 1.0) {
      Mm4R = Mm1R;
      Pm5R = Mm1R / MR;
    } else {
      Mm4R = Mm2R * (1.0 + 16.0 * beta * Mp2R);
      Pm5R = Mm2R * ((-2.0 - MR) + 16.0 * alpha * MR * Mp2R);
    }

    scrh = std::fmax(1.0 - sigma * Mm2, 0.0);
    M = Mp4L + Mm4R -
        Kp / fa * scrh * (pR - pL) / (rhoL + rhoR) / a / a *
            2.0; // also suspect
    m = a * M;
    m *= (M > 0.0) ? rhoL : rhoR;

    // Flux Calculation
    for (int var = 0; var < NumVar; ++var) {
      CellFlux[Tidx(var, i)] = (m > 0.0) ? m * PhiL[var] : m * PhiR[var];
    }
    CellFlux[Tidx(MOMX, i)] +=
        Pp5L * pL + Pm5R * pR -
        Ku * Pp5L * Pm5R * (rhoL + rhoR) * fa * a * (uR - uL);
  }
};
