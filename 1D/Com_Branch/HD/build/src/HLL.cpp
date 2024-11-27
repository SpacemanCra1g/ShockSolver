#include "../include/DomainClass.hpp"

void Domain::Hll(int Start, int Stop) {

  double SL, SR;
  double Lambda_Left_Minus, Lambda_Left_Plus;
  double Lambda_Right_Minus, Lambda_Right_Plus;
  double *LeftState_Prims, *RightState_Prims;
  int LeftIdx, RightIdx;

  Find_Cs(FluxWalls_Prims[LEFT], RS_CsL, Start, Stop);
  Find_Cs(FluxWalls_Prims[RIGHT], RS_CsR, Start, Stop);

  for (int i = Start; i < Stop; ++i) {

    LeftState_Prims = FluxWalls_Prims[RIGHT];
    RightState_Prims = FluxWalls_Prims[LEFT];

    LeftIdx = i;
    RightIdx = i + 1;

    SignalSpeed(LeftState_Prims, RS_CsR, LeftIdx, Lambda_Left_Minus,
                Lambda_Left_Plus);
    SignalSpeed(RightState_Prims, RS_CsL, RightIdx, Lambda_Right_Minus,
                Lambda_Right_Plus);

    SL = std::fmin(Lambda_Left_Minus, Lambda_Right_Minus);
    SR = std::fmax(Lambda_Left_Plus, Lambda_Right_Plus);

    // Left side Flux
    if (0.0 <= SL) {
      Flux(CellFlux, LeftState_Prims, LeftIdx, i);
    }
    // HLL Flux
    else if (0.0 <= SR) {
      HLL_Flux(CellFlux, LeftState_Prims, RightState_Prims, SL, SR, i);
    }
    // Right side Flux
    else {
      Flux(CellFlux, RightState_Prims, RightIdx, i);
    }
  }
};
