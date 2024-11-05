#include "../include/DomainClass.hpp"

double SIGN(double x) { return (x >= 0.0) ? 1.0 : -1.0; }

void Domain::Hllc(int Start, int Stop) {
  double SL, SR, Lam_CR, Lam_CL, Lam_RR, Lam_RL;
  double AL, BL, AR, BR, a, b, c, scrh, lamStar;
  double FL[5], FR[5], UL[5], UR[5], UStar[5];
  double p, div, dif;

  Cons2Prim(FluxWalls_Cons[LEFT], FluxWalls_Prims[LEFT], Start, Stop);
  Cons2Prim(FluxWalls_Cons[RIGHT], FluxWalls_Prims[RIGHT], Start, Stop);

  Find_Cs(FluxWalls_Prims[LEFT], RS_CsL, Start, Stop);
  Find_Cs(FluxWalls_Prims[RIGHT], RS_CsR, Start, Stop);
  for (int i = Start; i < Stop; ++i) {
    SignalSpeed(FluxWalls_Prims[LEFT], RS_CsL, i + 1, Lam_RL, Lam_RR);
    SignalSpeed(FluxWalls_Prims[RIGHT], RS_CsR, i, Lam_CL, Lam_CR);

    SL = std::fmin(Lam_CL, Lam_RL);
    SR = std::fmax(Lam_CR, Lam_RR);

    if (SL >= 0.0) {
      SR_Flux(CellFlux, FluxWalls_Prims[RIGHT], FluxWalls_Cons[RIGHT], i, i);
    } else if (SR <= 0.0) {
      SR_Flux(CellFlux, FluxWalls_Prims[LEFT], FluxWalls_Cons[LEFT], i + 1, i);
    } else {

      // Calculate Lam_Star
      SR_Flux(CellFlux, FluxWalls_Prims[RIGHT], FluxWalls_Cons[RIGHT], i, i);
      for (int var = 0; var < NumVar; ++var) {
        FL[var] = CellFlux[Tidx(var, i)];
        UL[var] = FluxWalls_Cons[RIGHT][Tidx(var, i)];
      }
      SR_Flux(CellFlux, FluxWalls_Prims[LEFT], FluxWalls_Cons[LEFT], i + 1, i);
      for (int var = 0; var < NumVar; ++var) {
        FR[var] = CellFlux[Tidx(var, i)];
        UR[var] = FluxWalls_Cons[LEFT][Tidx(var, i + 1)];
      }

      AL = SL * UL[ENER] - UL[MOMX];
      AR = SR * UR[ENER] - UR[MOMX];

      BL = UL[MOMX] * (SL - FluxWalls_Prims[RIGHT][Tidx(VELX, i)]) -
           FluxWalls_Prims[RIGHT][Tidx(PRES, i)];

      BR = UR[MOMX] * (SR - FluxWalls_Prims[LEFT][Tidx(VELX, i + 1)]) -
           FluxWalls_Prims[LEFT][Tidx(PRES, i + 1)];

      a = AR * SL - AL * SR;
      b = AL + BL * SR - AR - BR * SL;
      c = BR - BL;

      scrh = -0.5 * (b + SIGN(b) * sqrt(b * b - 4.0 * a * c));
      lamStar = c / scrh;

      p = (AL * lamStar - BL) / (1.0 - lamStar * SL);

      if (0.0 <= lamStar) {

        div = 1.0 / (SL - lamStar);
        dif = SL - FluxWalls_Prims[RIGHT][Tidx(VELX, i)];

        UStar[0] = UL[DENS] * dif * div;
        UStar[1] =
            (UL[MOMX] * dif + p - FluxWalls_Prims[RIGHT][Tidx(PRES, i)]) * div;
        UStar[2] = UL[MOMY] * dif * div;
        UStar[3] = UL[MOMZ] * dif * div;
        UStar[4] = (UL[ENER] * dif + p * lamStar -
                    FluxWalls_Prims[RIGHT][Tidx(PRES, i)] *
                        FluxWalls_Prims[RIGHT][Tidx(VELX, i)]) *
                   div;
        for (int var = 0; var < NumVar; ++var) {
          CellFlux[Tidx(var, i)] = SL * (UStar[var] - UL[var]) + FL[var];
        }
      } else {

        div = 1.0 / (SR - lamStar);
        dif = SR - FluxWalls_Prims[LEFT][Tidx(VELX, i + 1)];

        UStar[0] = UR[DENS] * dif * div;
        UStar[1] =
            (UR[MOMX] * dif + p - FluxWalls_Prims[LEFT][Tidx(PRES, i + 1)]) *
            div;
        UStar[2] = UR[MOMY] * dif * div;
        UStar[3] = UR[MOMZ] * dif * div;
        UStar[4] = (UR[ENER] * dif + p * lamStar -
                    FluxWalls_Prims[LEFT][Tidx(PRES, i + 1)] *
                        FluxWalls_Prims[LEFT][Tidx(VELX, i + 1)]) *
                   div;
        for (int var = 0; var < NumVar; ++var) {
          CellFlux[Tidx(var, i)] = SR * (UStar[var] - UR[var]) + FR[var];
        }
      }
    }
  }
}
