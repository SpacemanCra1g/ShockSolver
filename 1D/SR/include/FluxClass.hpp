#ifndef FLUXCLASS_H_
#define FLUXCLASS_H_

#include "../include/GP_Kernel.hpp"
#include "../include/Parameters.h"
#include <cmath>

class FluxClass {

public:
  GP_Kernel *Kern;
  double *Cons;
  double *dt;
  bool *Troubled;
  int *MoodOrd;
  bool MoodFinished = true;
  double *Uin;

  double FluxDir[NDIMS * 2][nqp][NumVar][xDim];

  double Flux[nqp][NumVar][xDim];

  void Fluxinit(double *COns, double *DT) {
    Cons = COns;
    dt = DT;
  }

  // Defined in the src/FOG.cpp file
  void FOG(int, int, int);

  // Defined in the src/GP-FVM.cpp file
  void GPR1(int, int, int);
  void GPR2(int, int, int);
  void GPR1Side(double *, int, int, int, int);
  void FOGSide(double *, int, int, int, int);
  void GPR2Side(double *, int, int, int, int);
  void Mood(int, int, int);

  // Defined in the src/WENO.cpp file
  void WENO(int, int, int);

  // Defined in the src/HLL.cpp file
  void HLL();
  void HLLSide(int);
  void SR_Flux(double *C, double *P, double *Flux);

  // Defined in the src/Detection.cpp file
  bool Detection(bool);

  void SpaceRecon() {
    for (int quad = 0; quad < nqp; ++quad) {
      for (int var = 0; var < NumVar; ++var) {
        for (int x = XStart - 2; x < XEnd + 2; ++x) {

#if SpaceMethod == Weno
          WENO(quad, var, x);
#elif SpaceMethod == Fog
          FOG(quad, var, x);
#elif SpaceMethod == Gp1
          GPR1(quad, var, x);
#elif SpaceMethod == Gp2
          GPR2(quad, var, x);
#elif SpaceMethod == Mood53
          Mood(quad, var, x);
#endif
        }
      }
    }
  }

  double GetPres(int x) {
    return (GAMMA - 1.0) * (Cons[Tidx(Ener, x)] -
                            0.5 * (Cons[Tidx(MomX, x)] * Cons[Tidx(MomX, x)] /
                                   Cons[Tidx(Dens, x)]));
  }

  double GetPresUin(int x) {
    return (GAMMA - 1.0) * (Uin[Tidx(Ener, x)] -
                            0.5 * (Uin[Tidx(MomX, x)] * Uin[Tidx(MomX, x)] /
                                   Uin[Tidx(Dens, x)]));
  }

  void Recon() {
    for (int quad = 0; quad < nqp; ++quad) {
      for (int var = 0; var < NumVar; ++var) {
        for (int x = XStart; x < XEnd; ++x) {
          Cons[Tidx(var, x)] -=
              *dt * (Flux[quad][var][x] - Flux[quad][var][x - 1]) / dx;
        }
      }
    }
  }

  void ReconSide(int quad, int var, int x) {
    Cons[Tidx(var, x)] =
        Uin[Tidx(var, x)] -
        *dt * (Flux[quad][var][x] - Flux[quad][var][x - 1]) / dx;
  }
};

#endif // FLUXCLASS_H_
