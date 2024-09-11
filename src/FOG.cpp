#include "../include/FluxClass.hpp"

void FluxClass::FOG() {

  for (int qp = 0; qp < nqp; ++qp) {
    for (int var = 0; var < NumVar; ++var) {
      for (int xdir = 0; xdir < xDim; ++xdir) {
        for (int ydir = 0; ydir < yDim; ++ydir) {
          if (xdir > 0) {
            FluxDir[Left][qp][var][idx(xdir, ydir)] =
                Cons[Tidx(var, xdir, ydir)];
          }

          if (xdir < xDim - 1) {
            FluxDir[Right][qp][var][idx(xdir, ydir)] =
                Cons[Tidx(var, xdir, ydir)];
          }

#if NDIMS > 1
          if (ydir < yDim - 1) {
            FluxDir[Bottom][qp][var][xdir * yDim + ydir] =
                Cons[Tidx(var, xdir, ydir)];
          }

          if (ydir > 0) {
            FluxDir[Top][qp][var][xdir * yDim + ydir] =
                Cons[Tidx(var, xdir, ydir)];
          }
#endif
        }
      }
    }
  }
}
