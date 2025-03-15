#include "../include/DomainClass.hpp"

void Domain::Hll(int Start, int Stop) {

  double SL, SR, scalar;
  double *LeftState_Prims, *RightState_Prims;
  double LFlux[NumVar], RFlux[NumVar];

  Find_Cs(FluxWalls_Prims[LEFT], RS_CsL, Start, Stop + 1);
  Find_Cs(FluxWalls_Prims[RIGHT], RS_CsR, Start, Stop);

  Prims2Cons(FluxWalls_Prims[LEFT], FluxWalls_Cons[LEFT], Start, Stop + 1);
  Prims2Cons(FluxWalls_Prims[RIGHT], FluxWalls_Cons[RIGHT], Start, Stop);

  for (int i = Start; i < Stop; ++i) {

    LeftState_Prims = FluxWalls_Prims[RIGHT];
    RightState_Prims = FluxWalls_Prims[LEFT];

    HLL_Speed(LeftState_Prims, RightState_Prims, RS_CsR, RS_CsL, i, SL, SR);

    FillFlux(LeftState_Prims, LFlux, RightState_Prims, RFlux, i);

    // Left side Flux
    if (0.0 <= SL) {
      for (int var = 0; var < NumVar; ++var) {
        CellFlux[Tidx(var, i)] = LeftState_Prims[Tidx(var, i)];
      }

    }
    // HLL Flux
    else if (0.0 <= SR) {

      for (int var = 0; var < NumVar; ++var) {
        scalar = SL * SR *
                 (FluxWalls_Cons[LEFT][Tidx(var, i + 1)] -
                  FluxWalls_Cons[RIGHT][Tidx(var, i)]);
        CellFlux[Tidx(var, i)] =
            (SR * LFlux[var] - SL * RFlux[var] + scalar) / (SR - SL);
      }

    }
    // Right side Flux
    else {
      for (int var = 0; var < NumVar; ++var) {
        CellFlux[Tidx(var, i)] = RightState_Prims[Tidx(var, i + 1)];
      }
    }
    // std::cout << CellFlux[Tidx(DENS, i)] << " " << Prims[Tidx(DENS, i)] << "
    // "
    //           << i << std::endl;
  }
};
