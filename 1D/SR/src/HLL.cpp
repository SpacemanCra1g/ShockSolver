#include "../include/DomainClass.hpp"
// #include "../include/SRVarConvert.hpp"

void Domain::Hll(int Start, int Stop) {

  double SL, SR;
  double Lambda_Left_Minus, Lambda_Left_Plus;
  double Lambda_Right_Minus, Lambda_Right_Plus;
  int TestVar = PRES;
  double *LeftState_Cons, *RightState_Cons;
  double *LeftState_Prims, *RightState_Prims;
  int LeftIdx, RightIdx;

  // int Start = XStart - 1;
  // int Stop = XEnd;

  Cons2Prim(FluxWalls_Cons[LEFT], FluxWalls_Prims[LEFT], Start, Stop);
  Cons2Prim(FluxWalls_Cons[RIGHT], FluxWalls_Prims[RIGHT], Start, Stop);

  Find_Cs(FluxWalls_Prims[LEFT], RS_CsL, Start, Stop);
  Find_Cs(FluxWalls_Prims[RIGHT], RS_CsR, Start, Stop);

  for (int i = Start; i < Stop; ++i) {
    std::cout << RS_CsL[i] << " " << i << std::endl;
  }

  std::cout << "\n \n" << std::endl;
  for (int i = Start; i < Stop; ++i) {
    std::cout << RS_CsR[i] << " " << i << std::endl;
  }
  // exit(0);

  for (int i = Start; i < Stop; ++i) {

    LeftState_Cons = FluxWalls_Cons[RIGHT];
    RightState_Cons = FluxWalls_Cons[LEFT];

    LeftState_Prims = FluxWalls_Prims[RIGHT];
    RightState_Prims = FluxWalls_Prims[LEFT];

    LeftIdx = i;
    RightIdx = i + 1;

    // Note the confusion between RS_CsR being associated with the left RP
    // The fluxwalls[RIGHT] and RS_CsR are Right sides of the CELLS, and
    // therefore the left side of the Riemann Problem. Confusing, but so it
    // goes.
    SignalSpeed(LeftState_Prims, RS_CsR, LeftIdx, Lambda_Left_Minus,
                Lambda_Left_Plus);
    SignalSpeed(RightState_Prims, RS_CsL, RightIdx, Lambda_Right_Minus,
                Lambda_Right_Plus);

    SL = std::fmin(Lambda_Left_Minus, Lambda_Right_Minus);
    SR = std::fmax(Lambda_Left_Plus, Lambda_Right_Plus);

    // std::cout << " Cell " << i << std::endl;
    // std::cout << "Left lambda " << RS_CsR[i] << " " << LambdaLL << " "
    //           << LambdaLR << " " << SL << std::endl;
    // std::cout << "Right lambda " << RS_CsL[i + 1] << " " << LambdaRL << " "
    //           << LambdaRR << " " << SR << std::endl;
    // std::cout << std::endl;

    // Left side Flux
    if (0.0 <= SL) {
      // std::cout << "LEFT CALLED";
      // exit(0);
      SR_Flux(CellFlux, LeftState_Prims, LeftState_Cons, LeftIdx, i);
    }
    // HLL Flux
    else if (0.0 <= SR) {
      SR_HLL_Flux(CellFlux, LeftState_Prims, LeftState_Cons, RightState_Prims,
                  RightState_Cons, SL, SR, i);
    }
    // Right side Flux
    else {
      // std::cout << "RIGHT CALLED";
      // exit(0);
      SR_Flux(CellFlux, RightState_Prims, RightState_Cons, RightIdx, i);
    }
  }
  // exit(0);
};

// void FluxClass::HLLSide(int xdir) {

//   int i, iPlus1, Dim, quadpoint;

//   double VelL, VelR, PresL, PresR, CsL, CsR, SL, SR, LF1, LF2, LF3, RF1, RF2,
//       RF3;

//   i = xdir;
//   iPlus1 = i + 1;
//   Dim = 0;
//   quadpoint = 0;

//   VelL = FluxDir[Dim + 1][quadpoint][Dim + 1][i] /
//          FluxDir[Dim + 1][quadpoint][Dens][i];

//   VelR = FluxDir[Dim][quadpoint][Dim + 1][iPlus1] /
//          FluxDir[Dim][quadpoint][Dens][iPlus1];

//   // This Pressure assumes 1D, make sure to change later
//   PresL = (GAMMA - 1.0) * (FluxDir[Dim + 1][quadpoint][Ener][i] -
//                            VelL * FluxDir[Dim + 1][quadpoint][MomX][i]
//                            / 2.0);

//   PresR = (GAMMA - 1.0) * (FluxDir[Dim][quadpoint][Ener][iPlus1] -
//                            VelR * FluxDir[Dim][quadpoint][MomX][iPlus1]
//                            / 2.0);

//   CsL = std::sqrt(GAMMA * PresL / FluxDir[Dim + 1][quadpoint][Dens][i]);
//   CsR = std::sqrt(GAMMA * PresR / FluxDir[Dim][quadpoint][Dens][iPlus1]);

//   SL = std::fmin(VelL - CsL, VelR - CsR);
//   SR = std::fmax(VelL + CsL, VelR + CsR);

//   if (0.0 <= SL) {
//     Flux[quadpoint][Dens][i] = FluxDir[Dim + 1][quadpoint][MomX][i];

//     Flux[quadpoint][MomX][i] =
//         FluxDir[Dim + 1][quadpoint][MomX][i] * VelL + PresL;

//     Flux[quadpoint][Ener][i] =
//         VelL * (FluxDir[Dim + 1][quadpoint][Ener][i] + PresL);

//   } else if (0.0 <= SR) {

//     LF1 = FluxDir[Dim + 1][quadpoint][MomX][i];

//     LF2 = FluxDir[Dim + 1][quadpoint][MomX][i] * VelL + PresL;

//     LF3 = VelL * (FluxDir[Dim + 1][quadpoint][Ener][i] + PresL);

//     RF1 = FluxDir[Dim][quadpoint][MomX][iPlus1];

//     RF2 = FluxDir[Dim][quadpoint][MomX][iPlus1] * VelL + PresL;

//     RF3 = VelL * (FluxDir[Dim][quadpoint][Ener][iPlus1] + PresL);

//     Flux[quadpoint][Dens][i] = (SR * LF1 - SL * RF1 +
//                                 SL * SR *
//                                     (FluxDir[Dim][quadpoint][Dens][iPlus1] -
//                                      FluxDir[Dim + 1][quadpoint][Dens][i])) /
//                                (SR - SL);

//     Flux[quadpoint][MomX][i] = (SR * LF2 - SL * RF2 +
//                                 SL * SR *
//                                     (FluxDir[Dim][quadpoint][MomX][iPlus1] -
//                                      FluxDir[Dim + 1][quadpoint][MomX][i])) /
//                                (SR - SL);

//     Flux[quadpoint][Ener][i] = (SR * LF3 - SL * RF3 +
//                                 SL * SR *
//                                     (FluxDir[Dim][quadpoint][Ener][iPlus1] -
//                                      FluxDir[Dim + 1][quadpoint][Ener][i])) /
//                                (SR - SL);

//   } else {

//     Flux[quadpoint][Dens][i] = FluxDir[Dim][quadpoint][MomX][iPlus1];

//     Flux[quadpoint][MomX][i] =
//         FluxDir[Dim][quadpoint][MomX][iPlus1] * VelR + PresR;

//     Flux[quadpoint][Ener][i] =
//         VelR * (FluxDir[Dim][quadpoint][Ener][iPlus1] + PresR);
//   }
// };
