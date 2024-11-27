#include "../include/DomainClass.hpp"

void Domain::Fog(int start, int stop) {
  double val;

  for (int xdir = start; xdir < stop; ++xdir) {
    for (int var = 0; var < NumVar; ++var) {
      val = Prims[Tidx(var, xdir)];
      FluxWalls_Prims[LEFT][Tidx(var, xdir)] = val;
      FluxWalls_Prims[RIGHT][Tidx(var, xdir)] = val;
    }
  }
}
