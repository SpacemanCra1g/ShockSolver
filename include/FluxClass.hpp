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

  double FluxDir[NDIMS * 2][nqp][NumVar][xDim * yDim];

  double Flux[nqp][NumVar][xDim * yDim];

  void Fluxinit(double *COns, double *DT) {
    Cons = COns;
    dt = DT;
  }

  // Defined in the src/FOG.cpp file
  void FOG(int, int, int, int);

  // Defined in the src/GP-FVM.cpp file
  void GP(int, int, int, int);
  void GPR1(int, int, int, int);
  void GPR2(int, int, int, int);
  void GPR1Side(int, int, int, int, int);
  void FOGSide(int, int, int, int, int);
  void GPR2Side(int, int, int, int, int);
  void Mood(int, int, int, int);

  // Defined in the src/WENO.cpp file
  void WENO(int, int, int, int);

  // Defined in the src/HLL.cpp file
  void HLL();

  // Defined in the src/Detection.cpp file
  bool Detection(bool);

  void SpaceRecon() {
    for (int quad = 0; quad < nqp; ++quad) {
      for (int var = 0; var < NumVar; ++var) {
        for (int x = 2; x < REdgeX - 2; ++x) {
          for (int y = YStart; y < YEnd; ++y) {
#if SpaceMethod == Weno
            WENO(quad, var, x, y);
#elif SpaceMethod == Fog
            FOG(quad, var, x, y);
#elif SpaceMethod == Gp1
            GPR1(quad, var, x, y);
#elif SpaceMethod == Gp2
            GPR2(quad, var, x, y);
#elif SpaceMethod == Mood53
            Mood(quad, var, x, y);
#endif
          }
        }
      }
    }
  }

  double GetPres(int x, int y) {
    return (GAMMA - 1.0) *
           (Cons[Tidx(Ener, x, y)] -
            0.5 * (Cons[Tidx(MomX, x, y)] * Cons[Tidx(MomX, x, y)] /
                   Cons[Tidx(Dens, x, y)]));
  }

  double GetPresUin(int x, int y) {
    return (GAMMA - 1.0) *
           (Uin[Tidx(Ener, x, y)] -
            0.5 * (Uin[Tidx(MomX, x, y)] * Uin[Tidx(MomX, x, y)] /
                   Uin[Tidx(Dens, x, y)]));
  }

  void Recon() {
    for (int quad = 0; quad < nqp; ++quad) {
      for (int var = 0; var < NumVar; ++var) {
        for (int x = XStart; x < XEnd; ++x) {
          for (int y = YStart; y < YEnd; ++y) {
            Cons[Tidx(var, x, y)] -=
                *dt *
                (Flux[quad][var][idx(x, y)] - Flux[quad][var][idx(x - 1, y)]) /
                dx;
          }
        }
      }
    }
  }

  void ReconSide(int quad, int var, int x, int y) {
    Cons[Tidx(var, x, y)] -=
        *dt * (Flux[quad][var][idx(x, y)] - Flux[quad][var][idx(x - 1, y)]) /
        dx;
  }
};

#endif // FLUXCLASS_H_
