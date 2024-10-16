#include "../include/DomainClass.hpp"
#include <cblas.h>
#include <cmath>

void Domain::NeumannBC() {
  for (int i = 0; i < NGC; ++i) {
    for (int var = 0; var < NumVar; ++var) {
      Cons[Tidx(var, i)] = Cons[Tidx(var, XStart)];
      Cons[Tidx(var, REdgeX - i - 1)] = Cons[Tidx(var, XEnd - 1)];
    }
  }
}
void Domain::ShuOsherBC() {
  for (int i = 0; i < NGC; ++i) {
    for (int var = 0; var < NDIMS; ++var) {
      // Left edge of X direction
      std::fill(&Cons[Tidx(var, 1)], &Cons[Tidx(var, NGC)], Cons[Tidx(var, 0)]);
      // Right edge of X direction
      std::fill(&Cons[Tidx(var, XEnd)], &Cons[Tidx(var, REdgeX)],
                Cons[Tidx(var, REdgeX - 1)]);
    }
  }
}
