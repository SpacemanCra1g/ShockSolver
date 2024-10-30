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
  double d = 3.857143;
  double vx = 2.629369;
  double p = 10.333333333;
  double e = p / ((GAMMA - 1.0) * d);

  for (int i = 0; i < NGC; ++i) {
    Cons[Tidx(DENS, i)] = d;
    Cons[Tidx(MOMX, i)] = d * vx;
    Cons[Tidx(ENER, i)] = vx * vx * d * 0.5 + e * d;

    for (int var = 0; var < NumVar; ++var) {
      Cons[Tidx(var, XEnd + i)] = Cons[Tidx(var, XEnd - 1)];
    }
  }
}
