#include "../include/FluxClass.hpp"

void FluxClass::HLL() {

  int i, iPlus1;

  double VelL, VelR, PresL, PresR, CsL, CsR, SL, SR, LF1, LF2, LF3, RF1, RF2,
      RF3;

  for (int Dim = 0; Dim < NDIMS * 2; Dim += 2) {
    for (int quadpoint = 0; quadpoint < nqp; ++quadpoint) {
      for (int xdir = 1; xdir < xDim - 1; ++xdir) {

        i = xdir;
        iPlus1 = i + 1;

        VelL = FluxDir[Dim + 1][quadpoint][Dim + 1][i] /
               FluxDir[Dim + 1][quadpoint][Dens][i];

        VelR = FluxDir[Dim][quadpoint][Dim + 1][iPlus1] /
               FluxDir[Dim][quadpoint][Dens][iPlus1];

        // This Pressure assumes 1D, make sure to change later
        PresL =
            (GAMMA - 1.0) * (FluxDir[Dim + 1][quadpoint][Ener][i] -
                             VelL * FluxDir[Dim + 1][quadpoint][MomX][i] / 2.0);

        PresR = (GAMMA - 1.0) *
                (FluxDir[Dim][quadpoint][Ener][iPlus1] -
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

          Flux[quadpoint][Dens][i] =
              (SR * LF1 - SL * RF1 +
               SL * SR *
                   (FluxDir[Dim][quadpoint][Dens][iPlus1] -
                    FluxDir[Dim + 1][quadpoint][Dens][i])) /
              (SR - SL);

          Flux[quadpoint][MomX][i] =
              (SR * LF2 - SL * RF2 +
               SL * SR *
                   (FluxDir[Dim][quadpoint][MomX][iPlus1] -
                    FluxDir[Dim + 1][quadpoint][MomX][i])) /
              (SR - SL);

          Flux[quadpoint][Ener][i] =
              (SR * LF3 - SL * RF3 +
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
