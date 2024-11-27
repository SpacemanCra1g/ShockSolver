#include "../include/DomainClass.hpp"
#include <cblas.h>

void Domain::ShockTubeIC() {
  for (int i = 0; i < REdgeX; ++i) {
    if ((i - NGC) * dx + dx * 0.5 < 0.5) {
      DensP[i] = 1.0;
      Pres[i] = 1.0;
      Xvel[i] = 0.0;
    } else {
      DensP[i] = 0.125;
      Pres[i] = 0.1;
      Xvel[i] = 0.0;
    }
  }
  Prims2Cons(Prims, Cons, 0, REdgeX);
}

void Domain::ShuOsherIC() {

  for (int i = 0; i < REdgeX; ++i) {
    if (dx * (i - XStart) <= 0.5) {
      DensP[i] = 3.857143;
      Xvel[i] = 2.629369;
      Pres[i] = 10.33333;
    } else {
      DensP[i] =
          1.0 + (0.2 * std::sin(5.0 * (-4.5 + (dx * 0.5 + dx * (i - XStart)))));
      Xvel[i] = 0.;
      Pres[i] = 1.;
    }
  }
  Prims2Cons(Prims, Cons, 0, REdgeX);
}
