#include "../include/FluxClass.hpp"

void FluxClass::GP(int qp, int var, int xdir, int ydir) {
#if SpaceMethod == Gp1
  GPR1(qp, var, xdir, ydir);

#elif SpaceMethod == Gp2
  GPR2(qp, var, xdir, ydir);
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

void FluxClass::GPR1Side(double *U, int qp, int var, int xdir, int ydir,
                         int face) {
  double valueLeft = 0.0;
  double valueRight = 0.0;
  for (int i = 0; i < 3; ++i) {

    valueLeft += U[Tidx(var, xdir - 1 + i, ydir)] * Kern->R1Left[i];
    valueRight += U[Tidx(var, xdir - 1 + i, ydir)] * Kern->R1Right[i];
  }

  // These are flipped from where they should be. Need to figure out
  if (face == Left) {
    FluxDir[Left][qp][var][idx(xdir, ydir)] = valueRight;
  } else if (face == Right) {
    FluxDir[Right][qp][var][idx(xdir, ydir)] = valueLeft;
  }
}

void FluxClass::GPR2Side(double *U, int qp, int var, int xdir, int ydir,
                         int face) {
  double valueLeft = 0.0;
  double valueRight = 0.0;

  for (int i = 0; i < 5; ++i) {
    valueLeft += U[Tidx(var, xdir - 2 + i, ydir)] * Kern->R2Left[i];
    valueRight += U[Tidx(var, xdir - 2 + i, ydir)] * Kern->R2Right[i];
  }

  //   These are flipped from where they should be. Need to figure out
  if (face == Right) {
    FluxDir[Right][qp][var][idx(xdir, ydir)] = valueLeft;
  } else if (face == Left) {
    FluxDir[Left][qp][var][idx(xdir, ydir)] = valueRight;
  } else {
    std::cout << "Failed ot set face";
    exit(0);
  }
}

void FluxClass::FOGSide(double *U, int qp, int var, int xdir, int ydir,
                        int face) {
  if (face == Left) {
    FluxDir[Right][qp][var][idx(xdir, ydir)] = U[Tidx(var, xdir, ydir)];
  } else if (face == Right) {
    FluxDir[Left][qp][var][idx(xdir, ydir)] = U[Tidx(var, xdir, ydir)];
  }
}
