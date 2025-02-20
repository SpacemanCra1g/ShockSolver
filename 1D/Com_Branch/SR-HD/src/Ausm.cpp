#include "../include/DomainClass.hpp"
void Domain::Ausm(int Start, int Stop) {

  double rL, uL, pL, aL, ML, HL, vL, wL;
  double rR, uR, pR, aR, MR, HR, vR, wR;
  double Mp, Pp;
  double Mm, Pm;
  double Fp[NumVar], Fm[NumVar];
  double PosF, NegF, lorR, lorL;

  // Find_Cs(FluxWalls_Prims[RIGHT], RS_CsL, Start, Stop);
  // Find_Cs(FluxWalls_Prims[LEFT], RS_CsR, Start, Stop + 1);

  for (int i = Start; i < Stop; ++i) {
    rL = FluxWalls_Prims[RIGHT][Tidx(DENSP, i)];
    uL = FluxWalls_Prims[RIGHT][Tidx(VELX, i)];
    vL = FluxWalls_Prims[RIGHT][Tidx(VELY, i)];
    wL = FluxWalls_Prims[RIGHT][Tidx(VELZ, i)];
    pL = FluxWalls_Prims[RIGHT][Tidx(PRES, i)];
    // pL = (pL < 0.0) ? 0.001 : pL;
    // rL = (rL < 0.0) ? 0.001 : rL;

    HL = 1.0 + (GAMMA / (GAMMA - 1.0)) * pL / rL;
    aL = std::sqrt(GAMMA * pL / (HL * rL));
    ML = uL / aL;

    rR = FluxWalls_Prims[LEFT][Tidx(DENSP, i + 1)];
    uR = FluxWalls_Prims[LEFT][Tidx(VELX, i + 1)];
    vR = FluxWalls_Prims[LEFT][Tidx(VELY, i + 1)];
    wR = FluxWalls_Prims[LEFT][Tidx(VELZ, i + 1)];
    pR = FluxWalls_Prims[LEFT][Tidx(PRES, i + 1)];
    // pR = (pR < 0.0) ? 0.001 : pR;
    // rR = (rR < 0.0) ? 0.001 : rR;

    HR = 1.0 + (GAMMA / (GAMMA - 1.0)) * pR / rR;
    aR = std::sqrt(GAMMA * pR / (HR * rR));
    MR = uR / aR;

    lorL = std::pow(1.0 - (uL * uL + vL * vL + wL * wL), -0.5);
    lorR = std::pow(1.0 - (uR * uR + vR * vR + wR * wR), -0.5);

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

    Fp[0] = PosF * lorL * rL;
    Fp[1] = PosF * rL * lorL * lorL * HL * uL + Pp;
    Fp[2] = PosF * rL * lorL * lorL * HL * vL;
    Fp[3] = PosF * rL * lorL * lorL * HL * wL;
    Fp[4] = PosF * rL * lorL * lorL * HL;

    Fm[0] = NegF * lorR * rR;
    Fm[1] = NegF * rR * lorR * lorR * HR * uR + Pm;
    Fm[2] = NegF * rR * lorR * lorR * HR * vR;
    Fm[3] = NegF * rR * lorR * lorR * HR * wR;
    Fm[4] = NegF * rR * lorR * lorR * HR;

    for (int var = 0; var < NumVar; ++var) {
      CellFlux[Tidx(var, i)] = Fp[var] + Fm[var];
    }
  }
};

void Domain::Autsm(int Start, int Stop) {

  int i;

  double aL, ML, MpL, PpL, asL2, asL, atL;
  double aR, MR, MmR, PmR, asR2, asR, atR;
  double a, m, mp, mm, p;
  double rhoL, pL, uL, rhoR, pR, uR;
  double HL, HR, Fp[NumVar], Fm[NumVar];
  double alpha = 3.0 / 16.0, beta = 0.125;
  double lorL, lorR, vL, wL, vR, wR;
  // double alpha = 0.0, beta = 0.0;

  for (i = Start; i < Stop; i++) {

    rhoL = FluxWalls_Prims[RIGHT][Tidx(DENS, i)];
    uL = FluxWalls_Prims[RIGHT][Tidx(VELX, i)];
    vL = FluxWalls_Prims[RIGHT][Tidx(VELY, i)];
    wL = FluxWalls_Prims[RIGHT][Tidx(VELZ, i)];
    pL = FluxWalls_Prims[RIGHT][Tidx(PRES, i)];
    pL = (pL < 0.0) ? 0.001 : pL;

    rhoR = FluxWalls_Prims[LEFT][Tidx(DENS, i + 1)];
    uR = FluxWalls_Prims[LEFT][Tidx(VELX, i + 1)];
    vR = FluxWalls_Prims[LEFT][Tidx(VELY, i + 1)];
    wR = FluxWalls_Prims[LEFT][Tidx(VELZ, i + 1)];
    pR = FluxWalls_Prims[LEFT][Tidx(PRES, i + 1)];
    pR = (pR < 0.0) ? 0.001 : pR;

    lorL = std::pow(1.0 - (uL * uL + vL * vL + wL * wL), -0.5);
    lorR = std::pow(1.0 - (uR * uR + vR * vR + wR * wR), -0.5);

    HL = 1.0 + (GAMMA / (GAMMA - 1.0)) * pL / rhoL;
    HR = 1.0 + (GAMMA / (GAMMA - 1.0)) * pR / rhoR;

    aL = std::sqrt(GAMMA * pL / (HL * rhoL));
    aR = std::sqrt(GAMMA * pR / (HR * rhoR));

    asL2 = uL * uL; // Works ok
    // asL2 = uL * uL + vL * vL + wL * wL;
    asL2 = aL * aL / (GAMMA - 1.0) + 0.5 * asL2;
    asL2 *= 2.0 * (GAMMA - 1.0) / (GAMMA + 1.0);

    asR2 = uR * uR;
    // asR2 = uR * uR + vR * vR + wR * wR;
    asR2 = aR * aR / (GAMMA - 1.0) + 0.5 * asR2;
    asR2 *= 2.0 * (GAMMA - 1.0) / (GAMMA + 1.0);

    asL = std::sqrt(asL2);
    asR = std::sqrt(asR2);

    atL = asL2 / std::fmax(asL, std::fabs(uL));
    atR = asR2 / std::fmax(asR, std::fabs(uR));

    a = std::fmin(atL, atR);
    a = std::sqrt(aL * aR);
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
    // HR = rhoR * (0.5 * uR * uR + pR / ((GAMMA - 1.0) * rhoR));
    // HL = rhoL * (0.5 * uL * uL + pL / ((GAMMA - 1.0) * rhoL));

    Fp[0] = mp * lorL * rhoL;
    Fp[1] = mp * rhoL * lorL * lorL * HL * uL + PpL * pL;
    Fp[2] = mp * rhoL * lorL * lorL * HL * vL;
    Fp[3] = mp * rhoL * lorL * lorL * HL * wL;
    Fp[4] = mp * rhoL * lorL * lorL * HL;

    Fm[0] = mm * lorR * rhoR;
    Fm[1] = mm * rhoR * lorR * lorR * HR * uR + PmR * pR;
    Fm[2] = mm * rhoR * lorR * lorR * HR * vR;
    Fm[3] = mm * rhoR * lorR * lorR * HR * wR;
    Fm[4] = mm * rhoR * lorR * lorR * HR;

    for (int var = 0; var < NumVar; ++var) {
      CellFlux[Tidx(var, i)] = Fp[var] + Fm[var];
    }
  }
};

void Domain::Autsmup(int Start, int Stop) {
  // This is now the Ausm+up solver

  int i;

  double aL, ML, MpL, PpL, asL2, asL, atL;
  double aR, MR, MmR, PmR, asR2, asR, atR;
  double a, m, mp, mm, p;
  double rhoL, pL, uL, rhoR, pR, uR;
  double HL, HR, Fp[3], Fm[3];
  double alpha = 3.0 / 16.0, beta = 0.125;
  double M2Bar, Mo2, fa;
  double Kp = .25, Ku = .75, sigma = 1.0, M2inf = .3;
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

    // atL = asL2 / std::fmax(asL, uL);
    // atR = asR2 / std::fmax(asR, -uR);

    a = std::fmin(atL, atR);
    /*
        a = 0.5*(aL + aR);
    */
    /* --------------------------------------------
            define split Mach numbers
            define pressure terms
       -------------------------------------------- */

    ML = uL / a;
    MR = uR / a;

    M2Bar = .5 * (uL * uL + uR * uR) / (a * a);
    Mo2 = std::fmin(1.0, std::fmax(M2Bar, M2inf));
    fa = std::sqrt(Mo2) * (2.0 - std::sqrt(Mo2));
    // fa = 1.0;

    // Define alpha
    alpha = (3.0 / 16.0) * (-4.0 + 5.0 * fa * fa);

    if (std::fabs(ML) >= 1.0) {
      MpL = 0.5 * (ML + std::fabs(ML));
      PpL = ML > 0.0 ? 1.0 : 0.0;
    } else {
      MpL = 0.25 * (ML + 1.0) * (ML + 1.0) +
            beta * (ML * ML - 1.0) * (ML * ML - 1.0);
      PpL = 0.25 * (ML + 1.0) * (ML + 1.0) * (2.0 - ML) +
            alpha * ML * (ML * ML - 1.0) * (ML * ML - 1.0);
    }

    if (std::fabs(MR) >= 1.0) {
      MmR = 0.5 * (MR - std::fabs(MR));
      PmR = MR > 0.0 ? 0.0 : 1.0;
    } else {
      MmR = -0.25 * (MR - 1.0) * (MR - 1.0) -
            beta * (MR * MR - 1.0) * (MR * MR - 1.0);
      PmR = 0.25 * (MR - 1.0) * (MR - 1.0) * (2.0 + MR) -
            alpha * MR * (MR * MR - 1.0) * (MR * MR - 1.0);
    }

    m = MpL + MmR -
        (Kp / fa) * std::fmax(1.0 - sigma * M2Bar, 0.0) * (pR - pL) /
            (a * a * (rhoL + rhoR) * .5);

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

    p = PpL * pL + PmR * pR -
        Ku * PpL * PmR * (rhoL + rhoR) * (fa * a) * (uR - uL);

    for (int var = 0; var < NumVar; ++var) {
      CellFlux[Tidx(var, i)] = Fp[var] + Fm[var];
    }
    CellFlux[Tidx(MOMX, i)] += p;
  }
};
