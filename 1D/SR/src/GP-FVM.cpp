#include "../include/DomainClass.hpp"

void Domain::GP1(int xdir) {
  double valueLeft;
  double valueRight;

  for (int var = 0; var < NumVar; ++var) {
    valueLeft = 0.0;
    valueRight = 0.0;

    for (int i = 0; i < 3; ++i) {
      valueLeft += Cons[Tidx(var, xdir - 1 + i)] * Ker.R1Left[i];
      valueRight += Cons[Tidx(var, xdir - 1 + i)] * Ker.R1Right[i];
    }

    // These are flipped from where they should be. Need to figure out
    FluxWalls_Cons[LEFT][Tidx(var, xdir)] = valueRight;
    FluxWalls_Cons[RIGHT][Tidx(var, xdir)] = valueLeft;
  }
}

void Domain::GP2(int xdir) {
  double valueLeft;
  double valueRight;

  for (int var = 0; var < NumVar; ++var) {
    valueLeft = 0.0;
    valueRight = 0.0;

    for (int i = 0; i < 5; ++i) {
      valueLeft += Cons[Tidx(var, xdir - 2 + i)] * Ker.R2Left[i];
      valueRight += Cons[Tidx(var, xdir - 2 + i)] * Ker.R2Right[i];
    }

    // These are flipped from where they should be. Need to figure out
    FluxWalls_Cons[LEFT][Tidx(var, xdir)] = valueRight;
    FluxWalls_Cons[RIGHT][Tidx(var, xdir)] = valueLeft;
  }
}

// void FluxClass::GPR1Side(double *U, int qp, int var, int xdir, int face) {
//   double valueLeft = 0.0;
//   double valueRight = 0.0;
//   for (int i = 0; i < 3; ++i) {

//     valueLeft += U[Tidx(var, xdir - 1 + i)] * Kern->R1Left[i];
//     valueRight += U[Tidx(var, xdir - 1 + i)] * Kern->R1Right[i];
//   }

//   // These are flipped from where they should be. Need to figure out
//   if (face == Left) {
//     FluxDir[Left][qp][var][xdir] = valueRight;
//   } else if (face == Right) {
//     FluxDir[Right][qp][var][xdir] = valueLeft;
//   }
// }

// void FluxClass::GPR2Side(double *U, int qp, int var, int xdir, int face) {
//   double valueLeft = 0.0;
//   double valueRight = 0.0;

//   for (int i = 0; i < 5; ++i) {
//     valueLeft += U[Tidx(var, xdir - 2 + i)] * Kern->R2Left[i];
//     valueRight += U[Tidx(var, xdir - 2 + i)] * Kern->R2Right[i];
//   }

//   //   These are flipped from where they should be. Need to figure out
//   if (face == Right) {
//     FluxDir[Right][qp][var][xdir] = valueLeft;
//   } else if (face == Left) {
//     FluxDir[Left][qp][var][xdir] = valueRight;
//   } else {
//     std::cout << "Failed ot set face";
//     exit(0);
//   }
// }

// void FluxClass::FOGSide(double *U, int qp, int var, int xdir, int face) {
//   if (face == Left) {
//     FluxDir[Right][qp][var][xdir] = U[Tidx(var, xdir)];
//   } else if (face == Right) {
//     FluxDir[Left][qp][var][xdir] = U[Tidx(var, xdir)];
//   }
// }
