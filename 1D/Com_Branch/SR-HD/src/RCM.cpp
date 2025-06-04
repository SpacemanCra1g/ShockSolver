#include "../include/DomainClass.hpp"
#include "../include/QuadExactSolver.hpp"
#include <random>
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

realkind Uniform(){
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(.0, 1.0);
    return (realkind)dis(gen);
}

void ValidState(realkind *State, double *Prims, int i){
  if (State[3] > 0.0) {return;}
  int Left = -1, Right = 1;
  while(Prims[Tidx(PRES,i+Left)] <= 0.0){
    Left--;
  }
  while(Prims[Tidx(PRES,i+Right)] <= 0.0){
    Right++;
  }
  double m = (double)(-Left) + (double)Right;
  double B = Prims[Tidx(PRES,i+Right)];
  double A = Prims[Tidx(PRES,i+Left)];
  double val = (B-A)/m;
  for(int p = 0; p <= (int)m; ++p){
    Prims[Tidx(PRES,i+Left+p)] = Prims[Tidx(PRES,i+Left)] + (double)p*val;
  }
  State[3] = (realkind)Prims[Tidx(PRES,i)];
}

void Domain::rcm(int Start, int Stop){
  
  realkind Seq;
  realkind StateL[4], StateR[4], Result[4];
  double d,vx,p,vy;
  realkind dxt = dx/dt;
  double alpha, lor, h;
  // std::cout << dx << std::endl;
  // exit(0);
  

  Seq = VanDerCorput(rcm_Counter); // .655
  // Seq = Uniform(); // .7407, .644, ,6755
  // Seq = (Seq > .5) ? Seq - 1.0 : Seq;
  Seq -= .5;

  #if RIEMANN != HYBRID
  rcm_Counter++;  
  #endif
  for (int i = Start+1; i < Stop; ++i){
    // rcm_Counter += 0;
    // Seq = VanDerCorput(rcm_Counter); // .655
    // Seq -= .5;
    // Seq = Uniform(); // .7407, .644, ,6755
    // Seq = (Seq > .5) ? Seq - 1.0 : Seq;
    // Seq *=.8;
    if (Seq > 0.0){                                             
      // Seq -= .5;
      StateL[0] = (realkind) FluxWalls_Prims[RIGHT][Tidx(DENSP, i-1)];
      StateR[0] = (realkind) FluxWalls_Prims[LEFT][Tidx(DENSP, i)];

      StateL[1] = (realkind) FluxWalls_Prims[RIGHT][Tidx(VELX, i-1)];
      StateR[1] = (realkind) FluxWalls_Prims[LEFT][Tidx(VELX, i)];

      StateL[2] = (realkind) FluxWalls_Prims[RIGHT][Tidx(VELY, i-1)];
      StateR[2] = (realkind) FluxWalls_Prims[LEFT][Tidx(VELY, i)];

      StateL[3] = (realkind) FluxWalls_Prims[RIGHT][Tidx(PRES, i-1)];
      StateR[3] = (realkind) FluxWalls_Prims[LEFT][Tidx(PRES, i)];
      ValidState(StateL, Prims, i-1);
      ValidState(StateR, Prims, i);
      
    }else{
      StateL[0] = (realkind) FluxWalls_Prims[RIGHT][Tidx(DENS, i)];
      StateR[0] = (realkind) FluxWalls_Prims[LEFT][Tidx(DENS, i+1)];

      StateL[1] = (realkind) FluxWalls_Prims[RIGHT][Tidx(VELX, i)];
      StateR[1] = (realkind) FluxWalls_Prims[LEFT][Tidx(VELX, i+1)];

      StateL[2] = (realkind) FluxWalls_Prims[RIGHT][Tidx(VELY, i)];
      StateR[2] = (realkind) FluxWalls_Prims[LEFT][Tidx(VELY, i+1)];

      StateL[3] = (realkind) FluxWalls_Prims[RIGHT][Tidx(PRES, i)];
      StateR[3] = (realkind) FluxWalls_Prims[LEFT][Tidx(PRES, i+1)];

      ValidState(StateL, Prims, i);
      ValidState(StateR, Prims, i+1);
    }
    // ExactSample(StateL, StateR, Result, Seq*dxt);
    // SolveRiemannFlux(StateL, StateR, Result, 0.0);

    // std::cout << "Cell Number = "  << i << " Seq = " << Seq*dxt << std::endl;

 
    double va = StateL[1]*StateL[1] + StateL[2]*StateL[2];
    if ( va > 1.0){
      // std::cout << "TRIGGERED 1";
      // exit(0);
      StateL[1] /= va;
      StateL[2] /= va;
    }

    va = StateR[1]*StateR[1] + StateR[2]*StateR[2];
    if ( va > 1.0){
      // std::cout << "TRIGGERED 2";
      // exit(0);
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