#include "../include/FluxClass.hpp"

void FluxClass::GP() {
#if SpaceMethod == Gp1
  for (int qp = 0; qp < nqp; ++qp) {
    for (int var = 0; var < NumVar; ++var) {
      for (int xdir = 1; xdir < xDim - 1; ++xdir) {
        for (int ydir = YStart; ydir < YEnd; ++ydir) {
          GPR1(qp, var, xdir, ydir);
        }
      }
    }
  }

#elif SpaceMethod == Gp2

  for (int qp = 0; qp < nqp; ++qp) {
    for (int var = 0; var < NumVar; ++var) {
      for (int xdir = 2; xdir < xDim - 2; ++xdir) {
        for (int ydir = YStart; ydir < YEnd; ++ydir) {
          GPR2(qp, var, xdir, ydir);
        }
      }
    }
  }

#endif
}

void FluxClass::GPR1(int qp, int var, int xdir, int ydir) {
  double valueLeft = 0.0;
  double valueRight = 0.0;
  for (int i = 0; i < 3; ++i) {
    valueLeft += Cons[Tidx(var, xdir - 1 + i, ydir)] * Kern->R1Left[i];
    valueRight += Cons[Tidx(var, xdir - 1 + i, ydir)] * Kern->R1Right[i];
  }

  // These are flipped from where they should be. Need to figure out
  FluxDir[Left][qp][var][idx(xdir, ydir)] = valueRight;
  FluxDir[Right][qp][var][idx(xdir, ydir)] = valueLeft;
}

void FluxClass::GPR2(int qp, int var, int xdir, int ydir) {
  double valueLeft = 0.0;
  double valueRight = 0.0;
  for (int i = 0; i < 5; ++i) {
    valueLeft += Cons[Tidx(var, xdir - 2 + i, ydir)] * Kern->R2Left[i];
    valueRight += Cons[Tidx(var, xdir - 2 + i, ydir)] * Kern->R2Right[i];
  }

  // These are flipped from where they should be. Need to figure out
  FluxDir[Right][qp][var][idx(xdir, ydir)] = valueLeft;
  FluxDir[Left][qp][var][idx(xdir, ydir)] = valueRight;
}
