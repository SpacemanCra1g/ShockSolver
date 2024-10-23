#include "../include/FluxClass.hpp"

void FluxClass::FOG(int qp, int var, int xdir) {
  if (xdir > 0) {
    FluxDir[Left][qp][var][xdir] = Cons[Tidx(var, xdir)];
  }

  if (xdir < xDim - 1) {
    FluxDir[Right][qp][var][xdir] = Cons[Tidx(var, xdir)];
  }
}
