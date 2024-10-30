#include "../include/DomainClass.hpp"

void Domain::Hll(int Start, int Stop) {

  double SL, SR;
  double Eig_CL, Eig_CR;
  double Eig_RL, Eig_RR;

  Cons2Prim(FluxWalls_Cons[LEFT], FluxWalls_Prims[LEFT], Start, Stop);
  Cons2Prim(FluxWalls_Cons[RIGHT], FluxWalls_Prims[RIGHT], Start, Stop);

  Find_Cs(FluxWalls_Prims[LEFT], RS_CsR, Start, Stop);
  Find_Cs(FluxWalls_Prims[RIGHT], RS_CsL, Start, Stop);

  for (int i = Start; i < Stop; ++i) {

    SignalSpeed(FluxWalls_Prims[LEFT], RS_CsR, i + 1, Eig_RL, Eig_RR);
    SignalSpeed(FluxWalls_Prims[RIGHT], RS_CsL, i, Eig_CL, Eig_CR);

    SL = std::fmin(Eig_RL, Eig_CL);
    SR = std::fmax(Eig_RR, Eig_CR);

    if (0.0 <= SL) {
      HD_Flux(CellFlux, FluxWalls_Prims[RIGHT], FluxWalls_Cons[RIGHT], i, i);
    }

    else if (0.0 <= SR) {
      HD_HLL_Flux(CellFlux, FluxWalls_Prims[RIGHT], FluxWalls_Cons[RIGHT],
                  FluxWalls_Prims[LEFT], FluxWalls_Cons[LEFT], SL, SR, i);
    }

    else {
      HD_Flux(CellFlux, FluxWalls_Prims[LEFT], FluxWalls_Cons[LEFT], i + 1, i);
    }
  }
};
