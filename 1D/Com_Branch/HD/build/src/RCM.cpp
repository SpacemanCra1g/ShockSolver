#include "../include/DomainClass.hpp"
#include "../include/HydroExact.hpp"
using namespace std;
double VanDerCorput(int i){
  double result = 0.0;
  double counter = 0.0;
  do{
    result += ((double)(i % 2 != 0))*std::pow(.5, counter + 1.0);
    counter++;
    i >>= 1;
  } while(i > 0);
  return result;
}

void Domain::rcm(int Start, int Stop){
  double Seq;
  double StateL[3], StateR[3], Result[3];
  double d,vx,p;
  double dxt = dx/dt;

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
    ExactSample(StateL, StateR, Result, Seq*dxt);

    d = Result[0];
    vx = Result[1];
    p = Result[2];

    Cons[Tidx(DENS, i)] = d;
    Cons[Tidx(MOMX, i)] = d * vx;
    Cons[Tidx(ENER, i)] = d * (vx * vx * 0.5 + p / ((GAMMA - 1.0) * d));
  }
}