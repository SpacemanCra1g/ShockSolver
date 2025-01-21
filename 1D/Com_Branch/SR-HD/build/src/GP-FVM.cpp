#include "../include/DomainClass.hpp"

void Domain::GP1(int start, int stop) {
  double valueLeft;
  double valueRight;

  for (int xdir = start; xdir < stop; ++xdir) {
    for (int var = 0; var < NumVar; ++var) {
      valueLeft = 0.0;
      valueRight = 0.0;

      for (int i = 0; i < 3; ++i) {
        valueLeft += Prims[Tidx(var, xdir - 1 + i)] * Ker.R1Left[i];
        valueRight += Prims[Tidx(var, xdir - 1 + i)] * Ker.R1Right[i];
      }

      // These are flipped from where they should be. Need to figure out
      FluxWalls_Prims[LEFT][Tidx(var, xdir)] = valueRight;
      FluxWalls_Prims[RIGHT][Tidx(var, xdir)] = valueLeft;
    }
  }
}

void Domain::GP2(int start, int stop) {
  double valueLeft;
  double valueRight;

  for (int xdir = start; xdir < stop; ++xdir) {
    for (int var = 0; var < NumVar; ++var) {
      valueLeft = 0.0;
      valueRight = 0.0;

      for (int i = 0; i < 5; ++i) {
        valueLeft += Prims[Tidx(var, xdir - 2 + i)] * Ker.R2Left[i];
        valueRight += Prims[Tidx(var, xdir - 2 + i)] * Ker.R2Right[i];
      }

      // These are flipped from where they should be. Need to figure out
      FluxWalls_Prims[LEFT][Tidx(var, xdir)] = valueRight;
      FluxWalls_Prims[RIGHT][Tidx(var, xdir)] = valueLeft;
    }
  }
}
