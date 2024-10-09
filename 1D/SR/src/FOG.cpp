#include "../include/FluxClass.hpp"

void FluxClass::FOG(int qp, int var, int xdir) {

  FluxDir[Left][qp][var][xdir] = Cons[Tidx(var, xdir)];

  FluxDir[Right][qp][var][xdir] = Cons[Tidx(var, xdir)];

#if NDIMS > 1
  if (ydir < yDim - 1) {
    FluxDir[Bottom][qp][var][xdir * yDim + ydir] = Cons[Tidx(var, xdir)];
  }

  if (ydir > 0) {
    FluxDir[Top][qp][var][xdir * yDim + ydir] = Cons[Tidx(var, xdir)];
  }
#endif
}
