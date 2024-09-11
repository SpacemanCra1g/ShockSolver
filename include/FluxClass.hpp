#ifndef FLUXCLASS_H_
#define FLUXCLASS_H_

#include "../include/Parameters.h"
#include <cmath>
//  This is maybe a mistake, but I'm going to define a Flux container to
//  hold all of flux data objects and methods

class FluxClass {

public:
  double *Cons;
  double *dt;

  double FluxDir[NDIMS * 2][nqp][NumVar][xDim * yDim];

  double Flux[nqp][NumVar][xDim * yDim];

  void Fluxinit(double *COns, double *DT) {
    Cons = COns;
    dt = DT;
  }

  // Defined in the src/FOG.cpp file
  void FOG();

  // Defined in the src/WENO.cpp file
  void WENO();

  // Defined in the src/HLL.cpp file
  void HLL();

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
};

#endif // FLUXCLASS_H_
