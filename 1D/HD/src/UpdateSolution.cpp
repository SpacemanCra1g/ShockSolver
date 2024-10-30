#include "../include/DomainClass.hpp"

void Domain::Recon(int start, int stop) {
  double Constant = dt / dx;
  for (int var = 0; var < NumVar; ++var) {
    for (int x = start; x < stop; ++x) {
      Cons[Tidx(var, x)] -=
          Constant * (CellFlux[Tidx(var, x)] - CellFlux[Tidx(var, x - 1)]);
    }
  }
}
