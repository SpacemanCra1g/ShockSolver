#include "../include/DomainClass.hpp"
#include <cblas.h>
#include <cmath>

void Domain::NeumannBC(std::string Vars) {

  if (Vars == "Cons") {
    for (int i = 0; i < NGC; ++i) {
      for (int var = 0; var < NumVar; ++var) {
        Cons[Tidx(var, i)] = Cons[Tidx(var, XStart)];
        Cons[Tidx(var, REdgeX - i - 1)] = Cons[Tidx(var, XEnd - 1)];
      }
    }
  } else if (Vars == "Prims") {
    for (int i = 0; i < NGC; ++i) {
      for (int var = 0; var < NumVar; ++var) {
        Prims[Tidx(var, i)] = Prims[Tidx(var, XStart)];
        Prims[Tidx(var, REdgeX - i - 1)] = Prims[Tidx(var, XEnd - 1)];
      }
    }
  }
}
void Domain::ShuOsherBC(std::string Vars) {

  double d = 3.857143;
  double vx = 2.629369;
  double p = 10.333333333;
  double e = p / ((GAMMA - 1.0) * d);

  for (int i = 0; i < NGC; ++i) {
    Cons[Tidx(Dens, i)] = d;
    Cons[Tidx(MomX, i)] = d * vx;
    Cons[Tidx(Ener, i)] = vx * vx * d * 0.5 + e * d;

    for (int var = 0; var < NumVar; ++var) {
      Cons[Tidx(var, XEnd + i)] = Cons[Tidx(var, XEnd - 1)];
    }
  }

  // if (Vars == "Cons") {
  //   for (int i = 0; i < NGC; ++i) {

  //     for (int var = 0; var < NDIMS; ++var) {
  //       // Left edge of X direction
  //       std::fill(&Cons[Tidx(var, 1)], &Cons[Tidx(var, NGC)],
  //                 Cons[Tidx(var, 0)]);

  //       // Right edge of X direction
  //       std::fill(&Cons[Tidx(var, XEnd)], &Cons[Tidx(var, REdgeX)],
  //                 Cons[Tidx(var, REdgeX - 1)]);
  //     }
  //   }
  // }
}
