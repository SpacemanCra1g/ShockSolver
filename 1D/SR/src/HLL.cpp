#include "../include/FluxClass.hpp"
#include "../include/GenVarConvert.hpp"

void FluxClass::HLL() {

  int i, iPlus1;

  double CsL, CsR, SL, SR;
  double PrimL[5];
  double PrimR[5];
  double ConL[5];
  double ConR[5];

  for (int Dim = 0; Dim < NDIMS * 2; Dim += 2) {
    for (int quadpoint = 0; quadpoint < nqp; ++quadpoint) {
      for (int xdir = 1; xdir < xDim - 1; ++xdir) {

        i = xdir;
        iPlus1 = i + 1;

        for (int var = 0; var < NumVar; ++var) {
          ConL[var] = FluxDir[Right][quadpoint][var][i];
          ConR[var] = FluxDir[Left][quadpoint][var][iPlus1];
        }

        Cons2Prim(ConL, PrimL);
        Cons2Prim(ConR, PrimR);

        CsL = SRHD_CS(ConL, PrimL);
        CsR = SRHD_CS(ConR, PrimR);

        // CsL = std::sqrt(GAMMA * PresL / FluxDir[Dim +
        // 1][quadpoint][Dens][i]); CsR = std::sqrt(GAMMA * PresR /
        // FluxDir[Dim][quadpoint][Dens][iPlus1]);p
        double LambdaLL, LambdaLR;
        double LambdaRL, LambdaRR;

        SignalSpeed(PrimL, CsL, LambdaLL, LambdaLR);
        SignalSpeed(PrimR, CsR, LambdaRL, LambdaRR);

        SL = std::fmin(LambdaLL, LambdaRL);
        SR = std::fmax(LambdaLR, LambdaRR);

        double FL[5];
        double FR[5];

        SR_Flux(ConL, PrimL, FL);
        SR_Flux(ConR, PrimR, FR);

        if (0.0 <= SL) {
          for (int var = 0; var < NumVar; ++var) {
            Flux[quadpoint][var][i] = FL[var];
          }

        } else if (0.0 <= SR) {
          for (int var = 0; var < NumVar; ++var) {
            Flux[quadpoint][var][i] = (SR * FL[var] - SL * FR[var] +
                                       SL * SR * (ConR[var] - ConL[var])) /
                                      (SR - SL);
          }

        } else {
          for (int var = 0; var < NumVar; ++var) {
            Flux[quadpoint][var][i] = FR[var];
          }
        }
      }
    }
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
