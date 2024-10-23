#include "../include/FluxClass.hpp"
#include "../include/VarConvert.hpp"
#include <variant>

void FluxClass::HLL() {

  double CsL, CsR, SL, SR, LamLR, LamLL, LamRL, LamRR;
  double PrimL[NumVar], PrimR[NumVar], ConL[NumVar], ConR[NumVar];
  double FL[NumVar], FR[NumVar];

  for (int Dim = 0; Dim < NDIMS * 2; Dim += 2) {
    int quadpoint = 0;
    for (int i = 1; i < xDim - 1; ++i) {

      for (int var = 0; var < NumVar; ++var) {
        ConL[var] = FluxDir[Right][quadpoint][var][i];
        ConR[var] = FluxDir[Left][quadpoint][var][i + 1];
      }

      ConConvert(ConL, PrimL);
      ConConvert(ConR, PrimR);

      CsL = HD_CS(PrimL);
      CsR = HD_CS(PrimR);

      SignalSpeed(PrimL, CsL, LamLL, LamLR);
      SignalSpeed(PrimR, CsR, LamRL, LamRR);

      SL = std::fmin(LamLL, LamRL);
      SR = std::fmax(LamLR, LamRR);

      if (0.0 <= SL) {
        FillFlux(PrimL, FL);

        for (int var = 0; var < NumVar; ++var) {
          Flux[quadpoint][var][i] = FL[var];
        }

      } else if (0.0 <= SR) {

        FillFlux(PrimL, FL);
        FillFlux(PrimR, FR);

        for (int var = 0; var < NumVar; ++var) {
          Flux[quadpoint][var][i] = (SR * FL[var] - SL * FR[var] +
                                     SL * SR * (ConR[var] - ConL[var])) /
                                    (SR - SL);
        }

      } else {

        FillFlux(PrimR, FR);

        for (int var = 0; var < NumVar; ++var) {
          Flux[quadpoint][var][i] = FR[var];
        }
      }
    }
  }
};

void FluxClass::HLLSide(int i) {

  double CsL, CsR, SL, SR, LamLR, LamLL, LamRL, LamRR;
  double PrimL[NumVar], PrimR[NumVar], ConL[NumVar], ConR[NumVar];
  double FL[NumVar], FR[NumVar];

  for (int Dim = 0; Dim < NDIMS * 2; Dim += 2) {
    int quadpoint = 0;

    for (int var = 0; var < NumVar; ++var) {
      ConL[var] = FluxDir[Right][quadpoint][var][i];
      ConR[var] = FluxDir[Left][quadpoint][var][i + 1];
    }

    ConConvert(ConL, PrimL);
    ConConvert(ConR, PrimR);

    CsL = HD_CS(PrimL);
    CsR = HD_CS(PrimR);

    SignalSpeed(PrimL, CsL, LamLL, LamLR);
    SignalSpeed(PrimR, CsR, LamRL, LamRR);

    SL = std::fmin(LamLL, LamRL);
    SR = std::fmax(LamLR, LamRR);

    if (0.0 <= SL) {
      FillFlux(PrimL, FL);

      for (int var = 0; var < NumVar; ++var) {
        Flux[quadpoint][var][i] = FL[var];
      }

    } else if (0.0 <= SR) {

      FillFlux(PrimL, FL);
      FillFlux(PrimR, FR);

      for (int var = 0; var < NumVar; ++var) {
        Flux[quadpoint][var][i] =
            (SR * FL[var] - SL * FR[var] + SL * SR * (ConR[var] - ConL[var])) /
            (SR - SL);
      }

    } else {

      FillFlux(PrimR, FR);

      for (int var = 0; var < NumVar; ++var) {
        Flux[quadpoint][var][i] = FR[var];
      }
    }
  }
};
