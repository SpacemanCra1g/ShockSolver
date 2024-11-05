#include "../include/DomainClass.hpp"

void Domain::Hll(int Start, int Stop) {

  double SL, SR;
  double Lambda_Left_Minus, Lambda_Left_Plus;
  double Lambda_Right_Minus, Lambda_Right_Plus;
  double *LeftState_Cons, *RightState_Cons;
  double *LeftState_Prims, *RightState_Prims;
  int LeftIdx, RightIdx;

  Cons2Prim(FluxWalls_Cons[LEFT], FluxWalls_Prims[LEFT], Start, Stop);
  Cons2Prim(FluxWalls_Cons[RIGHT], FluxWalls_Prims[RIGHT], Start, Stop);

  // while (err1) {
  //   for (int var = 0; var < NumVar; ++var) {
  //     FluxWalls_Cons[LEFT][Tidx(var, err1)] = Cons[Tidx(var, err1)];
  //   }
  //   err1 = Cons2Prim(FluxWalls_Cons[LEFT], FluxWalls_Prims[LEFT], err1,
  //   Stop);
  // }

  // while (err2) {
  //   for (int var = 0; var < NumVar; ++var) {
  //     FluxWalls_Cons[RIGHT][Tidx(var, err1)] = Cons[Tidx(var, err1)];
  //   }
  //   err2 = Cons2Prim(FluxWalls_Cons[LEFT], FluxWalls_Prims[LEFT], err2,
  //   Stop);
  // }

  // #if SpaceMethod == WENO
  //   bool Failure;
  //   for (int i = Start; i < Stop; ++i) {
  //     Failure = false;
  //     for (int var = 0; var < NumVar; ++var) {
  //       err1 = (FluxWalls_Cons[RIGHT][Tidx(var, i)] - Cons[Tidx(var, i)]) *
  //              (Cons[Tidx(var, i)] - FluxWalls_Cons[LEFT][Tidx(var, i)]);
  //       if (err1 < 0.0) {
  //         Failure = true;
  //       }
  //     }
  //     if (Failure) {
  //       for (int var = 0; var < NumVar; ++var) {
  //         FluxWalls_Cons[LEFT][Tidx(var, i)] = Cons[Tidx(var, i)];
  //       }
  //       err1 = Cons2Prim(FluxWalls_Cons[LEFT], FluxWalls_Prims[LEFT], i, i +
  //       1); err2 = Cons2Prim(FluxWalls_Cons[RIGHT], FluxWalls_Prims[RIGHT],
  //       i, i + 1);
  //     }
  //   }
  // #endif

  Find_Cs(FluxWalls_Prims[LEFT], RS_CsL, Start, Stop);
  Find_Cs(FluxWalls_Prims[RIGHT], RS_CsR, Start, Stop);

  for (int i = Start; i < Stop; ++i) {

    LeftState_Cons = FluxWalls_Cons[RIGHT];
    RightState_Cons = FluxWalls_Cons[LEFT];

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

      SR_Flux(CellFlux, LeftState_Prims, LeftState_Cons, LeftIdx, i);
    }
    // HLL Flux
    else if (0.0 <= SR) {
      SR_HLL_Flux(CellFlux, LeftState_Prims, LeftState_Cons, RightState_Prims,
                  RightState_Cons, SL, SR, i);
    }
    // Right side Flux
    else {

      SR_Flux(CellFlux, RightState_Prims, RightState_Cons, RightIdx, i);
    }
  }
};
