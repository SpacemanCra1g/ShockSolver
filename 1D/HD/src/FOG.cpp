#include "../include/DomainClass.hpp"

void Domain::Fog(int start, int stop) {
  for (int var = 0; var < NumVar; ++var) {
    for (int xdir = start; xdir < stop; ++xdir) {
      FluxWalls_Cons[LEFT][Tidx(var, xdir)] = Cons[Tidx(var, xdir)];
      FluxWalls_Cons[RIGHT][Tidx(var, xdir)] = Cons[Tidx(var, xdir)];
    }
  }
}
