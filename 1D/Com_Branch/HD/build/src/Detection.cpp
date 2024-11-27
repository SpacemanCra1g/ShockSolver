#include "../include/DomainClass.hpp"

bool Domain::Detection() {

  double Minlocal, Maxlocal;
  double UL, UR, UM, Mm;
  int count = 0;
  int x = XStart;

  bool DMP;

  double divV, gradP, pL, pR, threshold1, Ca, Mach;
  double D, d, mx, E, p, value;
  int LOrd, ROrd;

  MoodFinished = true;
  threshold1 = 5.0;

  std::fill(ConversionFailed, ConversionFailed + xDim, false);
  Cons2Prim(Cons, Prims, XStart - 1, XEnd + 1);
  Find_Cs(Prims, Cs, XStart - 1, XEnd + 1);

Detect:
  while (x < XEnd) {
    std::cout << "IN the While loop" << std::endl;
    exit(0);

    if (MoodOrd[x] == 1) {
      Troubled[x] = false;
      x++;
      goto Detect;
    }

    if (ConversionFailed[x]) {
      Troubled[x] = true;
      x++;
      goto Detect;
    }
    if (!Troubled[x]) {
      x++;
      goto Detect;
    }

    d = Prims[Tidx(DENSP, x)];

    p = Prims[Tidx(PRES, x)];
    D = Cons[Tidx(DENS, x)];
    mx = Cons[Tidx(MOMX, x)];

    E = Cons[Tidx(ENER, x)];

    if ((D <= 0.0) || (std::isnan(D))) {
      Troubled[x] = true;
      x++;
      goto Detect;
    }

    if (p <= 0.0 || std::isnan(p)) {
      Troubled[x] = true;
      x++;
      goto Detect;
    }

    // First shock detector
    divV = (Cons[Tidx(MOMX, x + 1)] - Cons[Tidx(MOMX, x - 1)]) / dx;
    divV *= 0.5 / Cons[Tidx(DENS, x)];
    divV -= 0.5 *
            (Cons[Tidx(MOMX, x)] *
             (Cons[Tidx(DENS, x + 1)] - Cons[Tidx(DENS, x - 1)]) / dx) /
            std::pow(Cons[Tidx(DENS, x)], 2);

    // Define Soundspeed
    Ca = std::sqrt(Cs[x]);

    Mach = std::pow(Cons[Tidx(MOMX, x)], 2) / pow(Cons[Tidx(DENS, x)], 2);
    Mach = std::sqrt(Mach / Ca);

    // Second Shock detector

    pL = Cons[Tidx(PRES, x - 1)];
    pR = Cons[Tidx(PRES, x + 1)];

    gradP = 0.5 * (std::fabs(pR - pL) / std::fmin(pL, pR)) / dx;
    if (divV < -dx * dx && Mach > 0.2 && gradP > threshold1) {
      // std::cout << "DPM Called" << std::endl;
      // exit(0);
      DMP = true;
    } else {
      DMP = false;
    }

    if (DMP) {
      std::cout << "DMP Called" << std::endl;
      UL = Cons[Tidx(DENS, x - 1)];
      UM = Cons[Tidx(DENS, x)];
      UR = Cons[Tidx(DENS, x + 1)];

      Minlocal = std::fmin(std::fmin(UL, UM), UR);
      Maxlocal = std::fmax(std::fmin(UL, UM), UR);

      Mm = Maxlocal - Minlocal;

      // Plateau detection
      if (Mm > dx * dx * dx) {
        if (Cons[Tidx(DENS, x)] > Maxlocal || Cons[Tidx(DENS, x)] < Minlocal) {
          Troubled[x] = true;
          x++;
          goto Detect;

          if (false) {
          }
        } else {
          Troubled[x] = false;
          x++;
          goto Detect;
        }

      } else {
        Troubled[x] = false;
        x++;
        goto Detect;
      }
    }

    for (int x = XStart - 1; x < XEnd + 1; ++x) {

      if (!Troubled[x]) {
        continue;
      }
      MoodFinished = false;

      MoodOrd[x] -= 2;
      if (MoodOrd[x] < 0) {
        std::cout << "1 called trouble" << std::endl;
        exit(0);
      }

      LOrd = std::fmin(MoodOrd[(idx(x - 1))], MoodOrd[x]);
      ROrd = std::fmin(MoodOrd[(idx(x + 1))], MoodOrd[x]);

      if (LOrd == 3) {
        for (int var = 0; var < NumVar; ++var) {
          value = 0.0;
          for (int j = 0; j < 3; ++j) {
            value += PrimsCopy[Tidx(var, x - 1 + j)] * Ker.R1Right[j];
          }
          FluxWalls_Prims[LEFT][Tidx(var, x)] = value;
        }
      } else if (LOrd == 1) {
        for (int var = 0; var < NumVar; ++var) {
          FluxWalls_Prims[LEFT][Tidx(var, x)] = PrimsCopy[Tidx(var, x)];
        }
      }
      if (ROrd == 3) {
        for (int var = 0; var < NumVar; ++var) {
          value = 0.0;
          for (int j = 0; j < 3; ++j) {
            value += PrimsCopy[Tidx(var, x - 1 + j)] * Ker.R1Left[j];
          }
          FluxWalls_Prims[RIGHT][Tidx(var, x)] = value;
        }
      } else if (ROrd == 1) {
        for (int var = 0; var < NumVar; ++var) {
          FluxWalls_Prims[RIGHT][Tidx(var, x)] = PrimsCopy[Tidx(var, x)];
        }
      }

      (this->*RiemannSolver)(x - 1, x + 2);
    }

    for (int x = XStart; x < XEnd; ++x) {

      if (!Troubled[x]) {
        continue;
      }

      Recon(x - 1, x + 1);
      ++count;
    }
  }
  std::cout << count << " Troubled Cells" << std::endl;
  return MoodFinished;
}
