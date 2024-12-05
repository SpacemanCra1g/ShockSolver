#include "../include/DomainClass.hpp"
#include <cmath>

// Note this Detection ONLY works for HD, Not SR_HD

bool Domain::U2Check(int x) {
  double *Center = &PrimsCopy[Tidx(DENS, x)];
  double ddx;

  // U2_MaxC[x] = -1.e14;
  // U2_MinC[x] = 1.e14;

  if (MoodOrd[x] == 5 || false) {
    for (int i = -2; i < 3; ++i) {
      Center = Center + i;
      ddx = *(Center - 1) + *(Center + 1) - 2 * (*Center);
      ddx /= dx * dx;
      U2_MaxC[x] = std::fmax(U2_MaxC[x], (ddx));
      U2_MinC[x] = std::fmin(U2_MinC[x], (ddx));
    }
  } else {
    for (int i = -1; i < 2; ++i) {
      Center = Center + i;
      ddx = *(Center - 1) + *(Center + 1) - 2 * (*Center);
      ddx /= dx * dx;
      U2_MaxC[x] = std::fmax(U2_MaxC[x], (ddx));
      U2_MinC[x] = std::fmin(U2_MinC[x], (ddx));
    }
  }
  return (U2_MaxC[x] * U2_MinC[x] > -dx &&
          std::fabs(U2_MinC[x] / U2_MaxC[x]) >= 0.5);
}

bool Domain::DMPCheck(int x) {
  double *Center, *MaxRho, *MinRho;

  Center = &PrimsCopy[Tidx(DENS, x)];
  MaxRho = &DMP_MaxRho[x];
  MinRho = &DMP_MinRho[x];

  *MaxRho = std::fmax(*MaxRho, *(Center + 1));
  *MaxRho = std::fmax(*MaxRho, *(Center - 1));

  *MinRho = std::fmin(*MinRho, *(Center + 1));
  *MinRho = std::fmin(*MinRho, *(Center - 1));

  return !(*MinRho <= Cons[Tidx(DENS, x)] && Cons[Tidx(DENS, x)] <= *MaxRho);
}

bool Domain::PlateauDetection(int x) {
  double MaxRho = -1.e14, MinRho = 1.e14;
  double *Center = &PrimsCopy[Tidx(DENS, x)];
  switch (MoodOrd[x]) {
  case 5:
    // MaxRho = std::fmax(*(Center - 2), *(Center - 1));
    // MaxRho = std::fmax(MaxRho, *(Center));
    // MaxRho = std::fmax(MaxRho, *(Center + 1));
    // MaxRho = std::fmax(MaxRho, *(Center + 2));

    // MinRho = std::fmin(*(Center - 2), *(Center - 1));
    // MinRho = std::fmin(MinRho, *(Center));
    // MinRho = std::fmin(MinRho, *(Center + 1));
    // MinRho = std::fmin(MinRho, *(Center + 2));
    for (int i = -2; i < 3; ++i) {
      MaxRho = std::fmax(*(Center + i), MaxRho);
      MinRho = std::fmin(*(Center + i), MinRho);
    }

    return (MaxRho - MinRho < dx * dx * dx);

  case 3:
    // MaxRho = std::fmax(*(Center), *(Center - 1));
    // MaxRho = std::fmax(MaxRho, *(Center + 1));

    // MinRho = std::fmin(*(Center), *(Center - 1));
    // MinRho = std::fmin(MinRho, *(Center + 1));

    for (int i = -1; i < 2; ++i) {
      MaxRho = std::fmax(*(Center + i), MaxRho);
      MinRho = std::fmin(*(Center + i), MinRho);
    }

    return (MaxRho - MinRho < dx * dx * dx);

  default:
    std::cout << "Error in Plateau Detection, MoodOrd not in bounds"
              << std::endl;
    exit(0);
  };
}

void Domain::ListofCellsToDetect(int start, int stop) {
  int x;
  IdxStop = 0;
  for (x = start; x < stop; ++x) {
    // On the first pass every cell is considered "troubled"
    // on subsequent passes, we only bother to run the detection on
    // previously troubled cells
    if (MoodOrd[x] == 1) {
      Troubled[x] = false;
      continue;
    }
    if (Troubled[x]) {
      if (Troubled[x - 1]) {
        Cons2Prim(Cons, Prims, x, x + 1);
        Find_Cs(Prims, Cs, x, x + 1);
      } else {
        Cons2Prim(Cons, Prims, x - 1, x + 1);
        Find_Cs(Prims, Cs, x - 1, x + 1);
      }
    } else {
      if (Troubled[x - 1]) {
        Cons2Prim(Cons, Prims, x, x + 1);
        Find_Cs(Prims, Cs, x, x + 1);
      } else {
        continue;
      }
    }
  }

  for (x = start; x < stop; ++x) {
    if (Troubled[x] || ConversionFailed[x]) {
      if (MoodOrd[x] <= 1) {
        Troubled[x] = false;
        continue;
      } else {
        // std::cout << "Stuck in this else lopp" << std::endl;
        // std::cout << "This is just getting goofy" << MoodOrd[x] << std::endl;
        // std::cout << "Cell Number " << x << std::endl;
        TroubledIdx[IdxStop] = x;
        IdxStop++;
      }
    }
  }
}

bool Domain::FirstShockDetector(int x) {
  double divV, gradP, pL, pR, Ca, Mach;
  double D, mx, p, d;

  D = Cons[Tidx(DENS, x)];
  mx = Cons[Tidx(MOMX, x)];
  d = Prims[Tidx(DENSP, x)];
  p = Prims[Tidx(PRES, x)];

  divV = (Cons[Tidx(MOMX, x + 1)] - Cons[Tidx(MOMX, x - 1)]) / dx;
  divV *= 0.5 / Cons[Tidx(DENS, x + 1)];
  divV -= 0.5 *
          (mx * (Cons[Tidx(DENS, x + 1)] - Cons[Tidx(DENS, x - 1)]) / dx) /
          (D * D);

  Ca = Cs[x];

  Mach = std::sqrt((mx * mx / (D * D)) / Ca);

  pL = Prims[Tidx(PRES, x - 1)];
  pR = Prims[Tidx(PRES, x + 1)];

  gradP = 0.5 * (std::fabs(pR - pL) / std::fmin(pL, pR) / dx);

  return (divV < -dx * dx && Mach > 0.2 && gradP > 5.0);
}

bool Domain::Detection() {

  int x, i;

  double D, d, p, value, threshold1, mx;

  int LOrd, ROrd;

  MoodFinished = false;
  threshold1 = 5.0;

  ListofCellsToDetect(XStart, XEnd);

  // for (i = 0; i < IdxStop; ++i) {
  //   std::cout << TroubledIdx[i] << std::endl;
  // }
  // exit(0);
  //
  // std::cout << "Numbe of Troubled Cells going into loop = " << IdxStop
  //           << std::endl;

  for (i = 0; i < IdxStop; ++i) {
    x = TroubledIdx[i];
    D = Cons[Tidx(DENS, x)];
    mx = Cons[Tidx(MOMX, x)];
    d = Prims[Tidx(DENSP, x)];
    p = Prims[Tidx(PRES, x)];

    // Positivity check
    if (d <= 0.0 || D <= 0.0 || p <= 0.0) {

      continue;
    }
    if (std::isnan(d) || std::isnan(p) || std::isnan(D)) {
      continue;
    }

    // We perform the CSD Checks in this step
    if (!FirstShockDetector(x)) {
      // std::cout << x << " Failed CSD Checkpoint" << std::endl;
      // exit(0);
      Troubled[x] = false;

      continue;
    }
    // If we are here then we have determined that we are inside a shock

    // Perform the Plateau Detection
    if (PlateauDetection(x)) {
      std::cout << x << " Failed Plateau Detection" << std::endl;
      Troubled[x] = false;

      continue;
    }
    // If we are here then we are not in a constant flat plateau state

    // Run DMP Check
    if (DMPCheck(x)) {
      // Run U2 Check, if fails then we know we have a truely troubled cell
      if (!U2Check(x) || true) {
        continue;
      }
    }
    // If either test failed the cell should not be cascaded
    Troubled[x] = false;

    continue;
  }

  // Reconstruction on our troubled cells

  IdxStop = 0;
  for (i = XStart; i < XEnd; ++i) {
    if (Troubled[i]) {
      if (MoodOrd[i] > 1) {
        TroubledIdx[IdxStop] = i;
        IdxStop++;
      }
    }
  }

  if (IdxStop == 0) {
    MoodFinished = true;
    return MoodFinished;
  } else {
    std::cout << IdxStop << " Troubled Cells on this pass" << std::endl;
  }

  // for (int p = 0; p < IdxStop; ++p) {
  //   std::cout << TroubledIdx[p] << std::endl;
  // }
  // exit(0);

  for (i = 0; i < IdxStop; i++) {
    x = TroubledIdx[i];

    MoodOrd[x] -= 2;
    if (MoodOrd[x] < 1) {
      // exit(0);
      MoodOrd[x] = 1;
      Troubled[x] = false;

      continue;
    }

    LOrd = std::fmin(MoodOrd[(idx(x - 1))], MoodOrd[x]);
    ROrd = std::fmin(MoodOrd[(idx(x + 1))], MoodOrd[x]);

    if (LOrd == 3) {
      Prims2Cons(PrimsCopy, Cons, x - 1, x + 2);
      for (int var = 0; var < NumVar; ++var) {
        value = 0.0;
        for (int j = 0; j < 3; ++j) {
          value += Cons[Tidx(var, x - 1 + j)] * Ker.R1Right[j];
        }
        FluxWalls_Cons[LEFT][Tidx(var, x)] = value;
      }
    } else if (LOrd == 1) {
      Prims2Cons(PrimsCopy, Cons, x, x + 1);
      for (int var = 0; var < NumVar; ++var) {
        FluxWalls_Cons[LEFT][Tidx(var, x)] = Cons[Tidx(var, x)];
      }
    }
    if (ROrd == 3) {
      Prims2Cons(PrimsCopy, Cons, x - 1, x + 2);
      for (int var = 0; var < NumVar; ++var) {
        value = 0.0;
        for (int j = 0; j < 3; ++j) {
          value += Cons[Tidx(var, x - 1 + j)] * Ker.R1Left[j];
        }
        FluxWalls_Cons[RIGHT][Tidx(var, x)] = value;
      }
    } else if (ROrd == 1) {
      for (int var = 0; var < NumVar; ++var) {
        FluxWalls_Cons[RIGHT][Tidx(var, x)] = Cons[Tidx(var, x)];
      }
    }

    Cons2Prim(FluxWalls_Cons[RIGHT], FluxWalls_Prims[RIGHT], x - 1, x + 2);
    Cons2Prim(FluxWalls_Cons[LEFT], FluxWalls_Prims[LEFT], x - 1, x + 2);

    (this->*RiemannSolver)(x - 1, x + 2);
    // Troubled[x - 1] = true;
    // Troubled[x] = true;
    // Troubled[x + 1] = true;
  }

  for (i = 0; i < IdxStop; i++) {
    x = TroubledIdx[i];
    Prims2Cons(PrimsCopy, Cons, x - 1, x + 2);
    Recon(x - 1, x + 2);
  }

  // std::cout << "Troubled Cells = " << IdxStop << " on this pass" <<
  // std::endl; std::cout << "Cell Number = " << TroubledIdx[IdxStop] <<
  // std::endl; std::cout << "Mood Ord = " << MoodOrd[TroubledIdx[IdxStop]] <<
  // std::endl; std::cout << "The time is " << T << std::endl;

  // if (IdxStop) {
  //   for (int var = 0; var < NumVar; ++var) {
  //     std::cout << Cons[Tidx(var, 4)] << " ";
  //   }
  //   exit(0);
  // }

  return MoodFinished;
}
