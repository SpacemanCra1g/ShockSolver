#include "../include/FluxClass.hpp"
#include "../include/SRVarConvert.hpp"

void FluxClass::HLL() {

  int i, iPlus1;

  double CsL, CsR, SL, SR;
  double PrimL[5];
  double PrimR[5];
  double ConL[5];
  double ConR[5];
  int quadpoint = 0;

  // for (int xdir = XStart - 1; xdir < XEnd; ++xdir) {
  //   for (int var = 0; var < NumVar; ++var) {
  //     std::cout << FluxDir[Left][quadpoint][var][xdir] << " ";
  //   }
  //   std::cout << xdir;
  //   std::cout << std::endl;
  // }
  // exit(0);

  for (int xdir = XStart - 1; xdir < XEnd; ++xdir) {

    i = xdir;
    iPlus1 = i + 1;

    for (int var = 0; var < NumVar; ++var) {
      ConL[var] = FluxDir[Right][quadpoint][var][i];
      ConR[var] = FluxDir[Left][quadpoint][var][iPlus1];
    }

    ConConvert(ConL, PrimL);
    ConConvert(ConR, PrimR);

    // for (int var = 0; var < NumVar; ++var) {
    //   std::cout << PrimR[var] << " ";
    // }
    // std::cout << i + 1;
    // std::cout << std::endl;

    CsL = SRH_CS(PrimL);
    CsR = SRH_CS(PrimR);

    // std::cout << CsL << " " << CsR << " " << i << std::endl;

    double LambdaLL, LambdaLR;
    double LambdaRL, LambdaRR;

    SignalSpeed(PrimL, CsL, LambdaLL, LambdaLR);
    SignalSpeed(PrimR, CsR, LambdaRL, LambdaRR);

    // std::cout << LambdaLL << " " << LambdaLR << " " << LambdaRL << " "
    //           << LambdaRR << " " << i;
    // std::cout << std::endl;

    SL = std::fmin(LambdaLL, LambdaRL);
    SR = std::fmax(LambdaLR, LambdaRR);

    // std::cout << SL << " " << SR << " " << i << std::endl;

    double FL[5];
    double FR[5];

    SR_Flux(ConL, PrimL, FL);
    SR_Flux(ConR, PrimR, FR);

    std::cout << "####################" << std::endl;
    std::cout << "Cell Number " << i << std::endl;

    std::cout << "Con Left" << "            ";
    for (int var = 0; var < NumVar; ++var) {
      std::cout << ConL[var] << " ";
    }
    std::cout << std::endl;

    std::cout << "Con Right" << "            ";
    for (int var = 0; var < NumVar; ++var) {
      std::cout << ConR[var] << " ";
    }
    std::cout << std::endl;

    std::cout << "Prim Left" << "            ";
    for (int var = 0; var < NumVar; ++var) {
      std::cout << PrimL[var] << " ";
    }
    std::cout << std::endl;

    std::cout << "Prim Right" << "            ";
    for (int var = 0; var < NumVar; ++var) {
      std::cout << PrimR[var] << " ";
    }
    std::cout << std::endl;

    std::cout << "Flux Left" << "            ";
    for (int var = 0; var < NumVar; ++var) {
      std::cout << FL[var] << " ";
    }
    std::cout << std::endl;

    std::cout << "Flux Right" << "            ";
    for (int var = 0; var < NumVar; ++var) {
      std::cout << FR[var] << " ";
    }

    std::cout << std::endl;
    std::cout << "SL" << "            ";
    std::cout << SL << std::endl;

    std::cout << "SR" << "            ";
    std::cout << SR << std::endl;

    if (0.0 <= SL) {
      for (int var = 0; var < NumVar; ++var) {
        Flux[quadpoint][var][i] = FL[var];
      }

    } else if (0.0 <= SR) {
      for (int var = 0; var < NumVar; ++var) {
        Flux[quadpoint][var][i] =
            (SR * FL[var] - SL * FR[var] + SL * SR * (ConR[var] - ConL[var])) /
            (SR - SL);
      }

    } else {
      for (int var = 0; var < NumVar; ++var) {
        Flux[quadpoint][var][i] = FR[var];
      }
    }

    std::cout << "FLUX " << "            ";

    for (int var = 0; var < NumVar; ++var) {
      std::cout << Flux[quadpoint][var][i] << " ";
    }
    std::cout << "\n\n\n" << std::endl;
  }
};

void FluxClass::HLLSide(int xdir) {

  int i, iPlus1, Dim, quadpoint;

  double VelL, VelR, PresL, PresR, CsL, CsR, SL, SR, LF1, LF2, LF3, RF1, RF2,
      RF3;

  i = xdir;
  iPlus1 = i + 1;
  Dim = 0;
  quadpoint = 0;

  VelL = FluxDir[Dim + 1][quadpoint][Dim + 1][i] /
         FluxDir[Dim + 1][quadpoint][Dens][i];

  VelR = FluxDir[Dim][quadpoint][Dim + 1][iPlus1] /
         FluxDir[Dim][quadpoint][Dens][iPlus1];

  // This Pressure assumes 1D, make sure to change later
  PresL = (GAMMA - 1.0) * (FluxDir[Dim + 1][quadpoint][Ener][i] -
                           VelL * FluxDir[Dim + 1][quadpoint][MomX][i] / 2.0);

  PresR = (GAMMA - 1.0) * (FluxDir[Dim][quadpoint][Ener][iPlus1] -
                           VelR * FluxDir[Dim][quadpoint][MomX][iPlus1] / 2.0);

  CsL = std::sqrt(GAMMA * PresL / FluxDir[Dim + 1][quadpoint][Dens][i]);
  CsR = std::sqrt(GAMMA * PresR / FluxDir[Dim][quadpoint][Dens][iPlus1]);

  SL = std::fmin(VelL - CsL, VelR - CsR);
  SR = std::fmax(VelL + CsL, VelR + CsR);

  if (0.0 <= SL) {
    Flux[quadpoint][Dens][i] = FluxDir[Dim + 1][quadpoint][MomX][i];

    Flux[quadpoint][MomX][i] =
        FluxDir[Dim + 1][quadpoint][MomX][i] * VelL + PresL;

    Flux[quadpoint][Ener][i] =
        VelL * (FluxDir[Dim + 1][quadpoint][Ener][i] + PresL);

  } else if (0.0 <= SR) {

    LF1 = FluxDir[Dim + 1][quadpoint][MomX][i];

    LF2 = FluxDir[Dim + 1][quadpoint][MomX][i] * VelL + PresL;

    LF3 = VelL * (FluxDir[Dim + 1][quadpoint][Ener][i] + PresL);

    RF1 = FluxDir[Dim][quadpoint][MomX][iPlus1];

    RF2 = FluxDir[Dim][quadpoint][MomX][iPlus1] * VelL + PresL;

    RF3 = VelL * (FluxDir[Dim][quadpoint][Ener][iPlus1] + PresL);

    Flux[quadpoint][Dens][i] = (SR * LF1 - SL * RF1 +
                                SL * SR *
                                    (FluxDir[Dim][quadpoint][Dens][iPlus1] -
                                     FluxDir[Dim + 1][quadpoint][Dens][i])) /
                               (SR - SL);

    Flux[quadpoint][MomX][i] = (SR * LF2 - SL * RF2 +
                                SL * SR *
                                    (FluxDir[Dim][quadpoint][MomX][iPlus1] -
                                     FluxDir[Dim + 1][quadpoint][MomX][i])) /
                               (SR - SL);

    Flux[quadpoint][Ener][i] = (SR * LF3 - SL * RF3 +
                                SL * SR *
                                    (FluxDir[Dim][quadpoint][Ener][iPlus1] -
                                     FluxDir[Dim + 1][quadpoint][Ener][i])) /
                               (SR - SL);

  } else {

    Flux[quadpoint][Dens][i] = FluxDir[Dim][quadpoint][MomX][iPlus1];

    Flux[quadpoint][MomX][i] =
        FluxDir[Dim][quadpoint][MomX][iPlus1] * VelR + PresR;

    Flux[quadpoint][Ener][i] =
        VelR * (FluxDir[Dim][quadpoint][Ener][iPlus1] + PresR);
  }
};

void FluxClass::SR_Flux(double *C, double *P, double *Flux) {
  Flux[Dens] = C[Dens] * P[VelX];
  Flux[MomX] = C[MomX] * P[VelX] + P[Pres];
  Flux[MomY] = C[MomY] * P[VelX];
  Flux[MomZ] = C[MomZ] * P[VelX];
  Flux[Ener] = C[MomX];
}
