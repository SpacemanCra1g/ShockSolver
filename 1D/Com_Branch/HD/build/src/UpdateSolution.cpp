#include "../include/DomainClass.hpp"

void Domain::Recon(int start, int stop) {
  for (int var = 0; var < NumVar; ++var) {
    for (int x = start; x < stop; ++x) {
      Cons[Tidx(var, x)] -=
          dt * (CellFlux[Tidx(var, x)] - CellFlux[Tidx(var, x - 1)]) / dx;
    }
  }
}
