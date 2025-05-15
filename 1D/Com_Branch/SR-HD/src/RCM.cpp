#include "../include/DomainClass.hpp"
#include "../include/QuadExactSolver.hpp"
using namespace std;
realkind VanDerCorput(int i){
  realkind result = 0.0;
  realkind counter = 0.0;
  do{
    result += (realkind)((double)(i % 2 != 0))*std::pow(.5, counter + 1.0);
    counter++;
    i >>= 1;
  } while(i > 0);
  return result;
}

void Domain::rcm(int Start, int Stop){
  
  realkind Seq;
  realkind StateL[4], StateR[4], Result[4];
  double d,vx,p,vy;
  realkind dxt = dx/dt;
  double alpha, lor, h;

  Seq = VanDerCorput(rcm_Counter);
  Seq = (Seq > .5) ? Seq - 1.0 : Seq;

  #if RIEMANN != HYBRID
  rcm_Counter++;  
  #endif

  for (int i = Start+1; i < Stop; ++i){
    if (Seq > 0.0){                                             
      StateL[0] = (realkind) FluxWalls_Prims[RIGHT][Tidx(DENSP, i-1)];
      StateR[0] = (realkind) FluxWalls_Prims[LEFT][Tidx(DENSP, i)];

      StateL[1] = (realkind) FluxWalls_Prims[RIGHT][Tidx(VELX, i-1)];
      StateR[1] = (realkind) FluxWalls_Prims[LEFT][Tidx(VELX, i)];

      StateL[2] = (realkind) FluxWalls_Prims[RIGHT][Tidx(VELY, i-1)];
      StateR[2] = (realkind) FluxWalls_Prims[LEFT][Tidx(VELY, i)];

      StateL[3] = (realkind) FluxWalls_Prims[RIGHT][Tidx(PRES, i-1)];
      StateR[3] = (realkind) FluxWalls_Prims[LEFT][Tidx(PRES, i)];
      
    }else{
      StateL[0] = (realkind) FluxWalls_Prims[RIGHT][Tidx(DENS, i)];
      StateR[0] = (realkind) FluxWalls_Prims[LEFT][Tidx(DENS, i+1)];

      StateL[1] = (realkind) FluxWalls_Prims[RIGHT][Tidx(VELX, i)];
      StateR[1] = (realkind) FluxWalls_Prims[LEFT][Tidx(VELX, i+1)];

      StateL[2] = (realkind) FluxWalls_Prims[RIGHT][Tidx(VELY, i)];
      StateR[2] = (realkind) FluxWalls_Prims[LEFT][Tidx(VELY, i+1)];

      StateL[3] = (realkind) FluxWalls_Prims[RIGHT][Tidx(PRES, i)];
      StateR[3] = (realkind) FluxWalls_Prims[LEFT][Tidx(PRES, i+1)];
    }
    // ExactSample(StateL, StateR, Result, Seq*dxt);
    // SolveRiemannFlux(StateL, StateR, Result, 0.0);

    // std::cout << "Cell Number = "  << i << " Seq = " << Seq*dxt << std::endl;
    double va = StateL[1]*StateL[1] + StateL[2]*StateL[2];
    if ( va > 1.0){
      StateL[1] /= va;
      StateL[2] /= va;
    }

    va = StateR[1]*StateR[1] + StateR[2]*StateR[2];
    if ( va > 1.0){
      StateR[1] /= va;
      StateR[2] /= va;
    }
    
    SolveRiemannFlux(StateL, StateR, Result, Seq*dxt);
    // std::cout << "Cell Number = "  << i << " Seq = " << Seq*dxt << std::endl;

    d = (double) Result[0];
    vx = (double) Result[1];
    vy = (double) Result[2];
    p = (double) Result[3];
    lor = 1.0/std::sqrt(1.0 - vy*vy - vx*vx);
    h = 1.0 + (GAMMA/(GAMMA-1.0))*p/d;
    alpha = d*h*lor*lor;

    Cons[Tidx(DENS, i)] = d*lor;
    Cons[Tidx(MOMX, i)] = alpha * vx;
    Cons[Tidx(MOMY, i)] = alpha * vy;
    Cons[Tidx(MOMZ, i)] = 0.0;
    Cons[Tidx(ENER, i)] = alpha - p;
  }
}