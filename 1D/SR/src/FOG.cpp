#include "../include/DomainClass.hpp"

void Domain::Fog(int xdir) {

  for (int var = 0; var < NumVar; ++var) {
    FluxWalls_Cons[LEFT][Tidx(var, xdir)] = Cons[Tidx(var, xdir)];
    FluxWalls_Cons[RIGHT][Tidx(var, xdir)] = Cons[Tidx(var, xdir)];
  }
}
