#include "../include/DomainClass.hpp"
#include "../include/QuadExactSolver.hpp"
using namespace std;
realkind VanDerCorput(int i){
  realkind result = 0.0;
  realkind counter = 0.0;
  do{
    result += ((double)(i % 2 != 0))*std::pow(.5, counter + 1.0);
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
  rcm_Counter++;  

  for (int i = Start+1; i < Stop; ++i){
    if (Seq > 0.0){
      for (int var = 0; var < 3; ++var){
      StateL[var] = FluxWalls_Prims[RIGHT][Tidx(var,i-1)];
      StateR[var] = FluxWalls_Prims[LEFT][Tidx(var,i)];
      }
    }else{
      for (int var = 0; var < 3; ++var){
        StateL[var] = FluxWalls_Prims[RIGHT][Tidx(var,i)];
        StateR[var] = FluxWalls_Prims[LEFT][Tidx(var,i+1)];
        }
    }
    // ExactSample(StateL, StateR, Result, Seq*dxt);
    SolveRiemannFlux(StateL, StateR, Result, 0.0);

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
    Cons[Tidx(ENER, i)] = alpha - p;
  }
}