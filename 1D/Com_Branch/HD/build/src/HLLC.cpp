#include "../include/DomainClass.hpp"
// Lots of credit to the PLUTO authors 

double SIGN(double x) { return (x >= 0.0) ? 1.0 : -1.0; }

void Domain::Hllc(int Start, int Stop) {
  double SL, SR, Lam_CR, Lam_CL, Lam_RR, Lam_RL;
  double FL[NumVar], FR[NumVar], UL[NumVar], UR[NumVar], UStar[NumVar];
  double qL, qR, wL, wR, vs, sch;
  double rhoL,rhoR,pL,pR,vL,vR, EL,ER,ML,MR, pStar, rhoLStar,rhoRStar;
  const double P = 1.3, D = 1.3;
  

  // Cons2Prim(FluxWalls_Cons[LEFT], FluxWalls_Prims[LEFT], Start, Stop);
  // Cons2Prim(FluxWalls_Cons[RIGHT], FluxWalls_Prims[RIGHT], Start, Stop);
  Prims2Cons(FluxWalls_Prims[LEFT], FluxWalls_Cons[LEFT], Start, Stop);
  Prims2Cons(FluxWalls_Prims[RIGHT], FluxWalls_Cons[RIGHT], Start, Stop);

  Find_Cs(FluxWalls_Prims[LEFT], RS_CsL, Start, Stop);
  Find_Cs(FluxWalls_Prims[RIGHT], RS_CsR, Start, Stop);
  for (int i = Start; i < Stop; ++i) {
    SignalSpeed(FluxWalls_Prims[LEFT], RS_CsL, i + 1, Lam_RL, Lam_RR);
    SignalSpeed(FluxWalls_Prims[RIGHT], RS_CsR, i, Lam_CL, Lam_CR);

    SL = std::fmin(Lam_CL, Lam_RL);
    SR = std::fmax(Lam_CR, Lam_RR);

    if (SL >= 0.0) {

      Flux(CellFlux, FluxWalls_Prims[RIGHT], i, i);

    } else if (SR <= 0.0) {
      Flux(CellFlux, FluxWalls_Prims[LEFT], i + 1, i);
    } else {

      rhoL = FluxWalls_Prims[RIGHT][Tidx(DENS,i)];
      vL = FluxWalls_Prims[RIGHT][Tidx(VELX,i)];
      pL = FluxWalls_Prims[RIGHT][Tidx(PRES,i)];

      rhoR = FluxWalls_Prims[LEFT][Tidx(DENS,i+1)];
      vR = FluxWalls_Prims[LEFT][Tidx(VELX,i+1)];
      pR = FluxWalls_Prims[LEFT][Tidx(PRES,i+1)];

      ML = FluxWalls_Cons[RIGHT][Tidx(MOMX,i)];
      EL = FluxWalls_Cons[RIGHT][Tidx(ENER,i)];

      MR = FluxWalls_Cons[LEFT][Tidx(MOMX,i+1)];
      ER = FluxWalls_Cons[LEFT][Tidx(ENER,i+1)];

      qL = -pL + ML*(SL - vL);
      qR = pR - MR*(SR - vR);

      wL = rhoL*(SL - vL);
      wR = -rhoR*(SR - vR);

      vs = (qR + qL)/(wR + wL);
      

      if (vs >= 0.0){

        sch = rhoL*(SL - vL)/(SL - vs);
        Flux(CellFlux, FluxWalls_Prims[RIGHT], i, i);
        UStar[0] = sch;
        UStar[1] = sch*vs;
        UStar[2] = sch*(EL/rhoL + (vs-vL)*(vs + pL/(rhoL*(SL-vL))));
        pStar = (UStar[2] / UStar[0] - UStar[1] * UStar[1] * 0.5) * (GAMMA - 1.0) * UStar[0];


        for (int var = 0; var < NumVar; ++var){
          CellFlux[Tidx(var,i)] += SL*(UStar[var] - FluxWalls_Cons[RIGHT][Tidx(var,i)]);
        }

      }else{
        
        sch = rhoR*(SR - vR)/(SR - vs);
        Flux(CellFlux, FluxWalls_Prims[LEFT], i+1, i);
        UStar[0] = sch;
        UStar[1] = sch*vs;
        UStar[2] = sch*(ER/rhoR + (vs-vR)*(vs + pR/(rhoR*(SL-vR))));
        pStar = (UStar[2] / UStar[0] - UStar[1] * UStar[1] * 0.5) * (GAMMA - 1.0) * UStar[0];
        for (int var = 0; var < NumVar; ++var){
          CellFlux[Tidx(var,i)] += SR*(UStar[var] - FluxWalls_Cons[LEFT][Tidx(var,i+1)]);
        }
        
      }
      // if (pStar / FluxWalls_Prims[RIGHT][Tidx(PRES,i)] > P){
      //   RcmReduction[i] = true;
      // }
      // if (pStar / FluxWalls_Prims[LEFT][Tidx(PRES,i+1)] > P){
      //    RcmReduction[i+1] = true;
      // }
      // rhoLStar =rhoL*(SL - vL)/(SL - vs);
      // rhoRStar =rhoR*(SR - vR)/(SR - vs);

      // VSSave[i] = std::fabs(rhoLStar/rhoRStar) - 1.0;

      // if (i > Start){
      //   if(VSSave[i-1] > D && vs > 0.0){
      //     RcmReduction[i] = true;
      //   }
      // }
      // if(std::fabs(rhoLStar/rhoRStar - 1.0) > D && vs < 0.0){
      //     RcmReduction[i] = true;
      // }
    }
  }
}
