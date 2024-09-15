#include "../include/FluxClass.hpp"

void FluxClass::FOG(int qp, int var, int xdir, int ydir) {

  if (xdir > 0) {
    FluxDir[Left][qp][var][idx(xdir, ydir)] = Cons[Tidx(var, xdir, ydir)];
  }

  if (xdir < xDim - 1) {
    FluxDir[Right][qp][var][idx(xdir, ydir)] = Cons[Tidx(var, xdir, ydir)];
  }

#if NDIMS > 1
  if (ydir < yDim - 1) {
    FluxDir[Bottom][qp][var][xdir * yDim + ydir] = Cons[Tidx(var, xdir, ydir)];
  }

  if (ydir > 0) {
    FluxDir[Top][qp][var][xdir * yDim + ydir] = Cons[Tidx(var, xdir, ydir)];
  }
#endif
}
