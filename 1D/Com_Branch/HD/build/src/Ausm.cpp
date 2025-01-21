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
    aL = std::sqrt(GAMMA * pL / rL);
    ML = uL / aL;
    HL = aL * aL / (GAMMA - 1.0) + 0.5 * uL * uL;
    // HL = (rL * (uL * uL * .5 + pL / ((GAMMA - 1.0) * rL)) + pL) / rL;

    rR = FluxWalls_Prims[LEFT][Tidx(DENSP, i + 1)];
    uR = FluxWalls_Prims[LEFT][Tidx(VELX, i + 1)];
    pR = FluxWalls_Prims[LEFT][Tidx(PRES, i + 1)];
    aR = std::sqrt(GAMMA * pR / rR);
    MR = uR / aR;
    HR = aR * aR / (GAMMA - 1.0) + 0.5 * uR * uR;
    // HR = (rR * (uR * uR * .5 + pR / ((GAMMA - 1.0) * rR)) + pR) / rR;

    if (ML <= -1.0) {
      Mp = 0.0;
      Pp = 0.0;
    } else if (ML < 1.0) {
      Mp = (ML + 1.0) * (ML + 1.0) * 0.25;
      Pp = pL * (1.0 + ML) * (1.0 + ML) * (2.0 - ML) * 0.25;
      // Pp = (1.0 + ML) * pL * .5;
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
      // Pm = (1.0 - MR) * pR * 0.5;
    } else {
      Mm = 0.0;
      Pm = 0.0;
    }

    PosF = std::fmax(0.0, Mp + Mm) * aL;
    NegF = std::fmin(0.0, Mp + Mm) * aR;

    // if (i == 202 || i == 203) {
    //   std::cout << "PosF = " << PosF << " NegF = " << NegF << std::endl;
    // }

    Fp[0] = PosF * rL;
    Fp[1] = PosF * rL * uL + Pp;
    Fp[2] = PosF * rL * HL;

    Fm[0] = NegF * rR;
    Fm[1] = NegF * rR * uR + Pm;
    Fm[2] = NegF * rR * HR;

    for (int var = 0; var < NumVar; ++var) {
      CellFlux[Tidx(var, i)] = Fp[var] + Fm[var];
    }

    // if (i == 202 || i == 203) {
    //   std::cout << "CellFlux = " << CellFlux[Tidx(VELX, i)] << std::endl;
    // }
  }
};

void Domain::Autsm(int Start, int Stop) {

  double rL, uL, pL, aL, ML, HL, aStrL;
  double rR, uR, pR, aR, MR, HR, aStrR;
  double Mp, Pp;
  double Mm, Pm;
  double Fp[3], Fm[3];
  double PosF, NegF;
  double A;
  double alpha = 0.0; // 3.0 / 16.0;
  double beta = 0.0;  // 1.0 / 8.0;

  for (int i = Start; i < Stop; ++i) {
    rL = FluxWalls_Prims[RIGHT][Tidx(DENSP, i)];
    uL = FluxWalls_Prims[RIGHT][Tidx(VELX, i)];
    pL = FluxWalls_Prims[RIGHT][Tidx(PRES, i)];
    aL = std::sqrt(GAMMA * pL / rL);
    HL = aL * aL / (GAMMA - 1.0) + 0.5 * uL * uL;
    aStrL = std::sqrt(2.0 * (GAMMA - 1.0) * HL / (GAMMA + 1.0));

    rR = FluxWalls_Prims[LEFT][Tidx(DENSP, i + 1)];
    uR = FluxWalls_Prims[LEFT][Tidx(VELX, i + 1)];
    pR = FluxWalls_Prims[LEFT][Tidx(PRES, i + 1)];
    aR = std::sqrt(GAMMA * pR / rR);
    HR = aR * aR / (GAMMA - 1.0) + 0.5 * uR * uR;
    aStrR = std::sqrt(2.0 * (GAMMA - 1.0) * HR / (GAMMA + 1.0));

    A = std::fmin(aStrL * aStrL / std::fmax(aStrL, std::fabs(uL)),
                  aStrR * aStrR / std::fmax(aStrR, std::fabs(uR)));

    ML = uL / A;
    MR = uR / A;

    if (std::fabs(ML) >= 1.0) {
      Mp = 0.5 * (ML + std::fabs(ML));
      Pp = .5 * (1.0 + ML / std::fabs(ML));
    } else {
      Mp = 0.5 * (ML + 1.0) * (ML + 1.0) +
           beta * (ML * ML - 1.0) * (ML * ML - 1.0);
      Pp = .25 * (ML + 1.0) * (ML + 1.0) * (2.0 - ML) +
           alpha * ML * (ML * ML - 1.0) * (ML * ML - 1.0);
    }

    if (std::fabs(MR) >= 1.0) {
      Mm = .5 *
           (MR - std::fabs(MR)); // Code calls this a +, contrary to the paper
      Pm = .5 * (1.0 - MR / std::fabs(MR));
    } else {
      Mm = -0.5 * (MR - 1.0) * (MR - 1.0) -
           beta * (MR * MR - 1.0) * (MR * MR - 1.0);
      Pm = .25 * (MR - 1.0) * (MR - 1.0) * (2.0 + MR) -
           alpha * MR * (MR * MR - 1.0) * (MR * MR - 1.0);
      // Pm = (1.0 - MR) * pR * 0.5;
    }

    PosF = .5 * (Mp + Mm + std::fabs(Mp + Mm)) * A;
    NegF = .5 * (Mp + Mm - std::fabs(Mp + Mm)) * A;

    // if (i == 202 || i == 203) {
    //   std::cout << "PosF = " << PosF << " NegF = " << NegF << std::endl;
    // }

    Fp[0] = PosF * rL;
    Fp[1] = PosF * rL * uL + Pp * pL;
    Fp[2] = PosF * HL * rL; // Code doesnt include the rho?

    Fm[0] = NegF * rR;
    Fm[1] = NegF * rR * uR + Pm * pR;
    Fm[2] = NegF * rR * HR; // Code doesnt include the rho?

    for (int var = 0; var < NumVar; ++var) {
      CellFlux[Tidx(var, i)] = Fp[var] + Fm[var];
    }

    // if (i == 202 || i == 203) {
    //   std::cout << "CellFlux = " << CellFlux[Tidx(VELX, i)] << std::endl;
    // }
  }
};
