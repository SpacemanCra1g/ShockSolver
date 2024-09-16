#include "../include/FluxClass.hpp"

bool FluxClass::Detection(bool MoodFinished) {

  double Minlocal, Maxlocal, p;
  double UL, UR, UM, Mm;
  int count = 0;

  bool DMP;

  int k, ord, ord_derr;

  double divV, gradP, pL, pR, pT, pB, threshold1, threshold2, Ca, Mach;

  MoodFinished = true;

  // initialize DMP here

  threshold1 = 5.0;

  threshold2 = threshold1;

  for (int x = 2; x < REdgeX - 2; ++x) {
    for (int y = YStart; y < YEnd; ++y) {
      if (MoodOrd[idx(x, y)] == 1) {
        continue;
      }

      if ((Cons[Tidx(Dens, x, y)] <= 0.) ||
          (std::isnan(Cons[Tidx(Dens, x, y)]))) {
        Troubled[idx(x, y)] = true;
        continue;
      }

      p = GetPres(x, y);

      if (p <= 0.0 || std::isnan(p)) {
        Troubled[idx(x, y)] = true;
        continue;
      }

      // First shock detector
      divV = (Uin[Tidx(MomX, x + 1, y)] - Uin[Tidx(MomX, x - 1, y)]) / dx;
      divV *= 0.5 / Uin[Tidx(Dens, x, y)];
      divV -= 0.5 *
              (Uin[Tidx(MomX, x, y)] *
               (Uin[Tidx(Dens, x + 1, y)] - Uin[Tidx(Dens, x - 1, y)]) / dx) /
              std::pow(Uin[Tidx(Dens, x, y)], 2);

      // Define Soundspeed
      Ca = GetPresUin(x, y) * GAMMA / Uin[(Tidx(Dens, x, y))];

      Mach = std::pow(Uin[Tidx(MomX, x, y)], 2) / pow(Uin[Tidx(Dens, x, y)], 2);
      Mach = std::sqrt(Mach / Ca);

      // Second Shock detector

      pL = GetPresUin(x - 1, y);
      pR = GetPresUin(x + 1, y);

      gradP = 0.5 * (std::fabs(pR - pL) / std::fmin(pL, pR)) / dx;
      if (divV < -dx * dx && Mach > 0.2 && gradP > threshold1) {
        // std::cout << "DPM Called" << std::endl;
        // exit(0);
        DMP = true;
      } else {
        DMP = false;
      }

      if (DMP) {
        UL = Uin[Tidx(Dens, x - 1, y)];
        UM = Uin[Tidx(Dens, x, y)];
        UR = Uin[Tidx(Dens, x + 1, y)];

        Minlocal = std::fmin(std::fmin(UL, UM), UR);
        Maxlocal = std::fmax(std::fmin(UL, UM), UR);

        Mm = Maxlocal - Minlocal;

        // Plateau detection
        if (Mm > dx * dx * dx) {
          if (Cons[Tidx(Dens, x, y)] > Maxlocal ||
              Cons[Tidx(Dens, x, y)] < Minlocal) {
            Troubled[idx(x, y)] = true;

            if (false) {
            }
          }
        }
      }
    }
  }

  for (int x = 2; x < REdgeX - 2; ++x) {
    for (int y = YStart; y < YEnd; ++y) {
      if (!Troubled[idx(x, y)]) {
        continue;
      }
      MoodFinished = false;

      MoodOrd[(idx(x, y))] -= 2;
      if (MoodOrd[(idx(x, y))] < 0) {
        std::cout << "1 called trouble" << std::endl;
        exit(0);
      }

      int LOrd, ROrd;
      LOrd = std::fmin(MoodOrd[(idx(x - 1, y))], MoodOrd[(idx(x, y))]);
      ROrd = std::fmin(MoodOrd[(idx(x + 1, y))], MoodOrd[(idx(x, y))]);

      if (LOrd == 3) {
        for (int var = 0; var < NumVar; ++var) {
          GPR1Side(Uin, 0, var, x, y, Left);
          GPR1Side(Uin, 0, var, x - 1, y, Right);
          HLLSide(x - 1, y);
        }
      } else if (LOrd == 1) {
        for (int var = 0; var < NumVar; ++var) {
          FOGSide(Uin, 0, var, x, y, Left);
          FOGSide(Uin, 0, var, x - 1, y, Right);
          HLLSide(x - 1, y);
        }
      }
      if (ROrd == 3) {
        for (int var = 0; var < NumVar; ++var) {
          GPR1Side(Uin, 0, var, x, y, Right);
          GPR1Side(Uin, 0, var, x + 1, y, Left);
          HLLSide(x + 1, y);
        }
      } else if (ROrd == 1) {
        for (int var = 0; var < NumVar; ++var) {
          FOGSide(Uin, 0, var, x, y, Right);
          FOGSide(Uin, 0, var, x + 1, y, Left);
          HLLSide(x + 1, y);
        }
      }
    }
  }

  for (int x = 2; x < REdgeX - 2; ++x) {
    for (int y = YStart; y < YEnd; ++y) {
      if (!Troubled[idx(x, y)]) {
        continue;
      }
      for (int var = 0; var < NumVar; ++var) {
        ReconSide(0, var, x, y);
        ReconSide(0, var, x - 1, y);
        ReconSide(0, var, x + 1, y);
      }
      Troubled[idx(x, y)] = false;
      ++count;
    }
  }
  std::cout << count << " Troubled Cells" << std::endl;
  return MoodFinished;
};
