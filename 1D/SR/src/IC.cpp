#include "../include/DomainClass.hpp"
#include <cblas.h>

void Domain::ShockTubeIC() {
  for (int i = 0; i < REdgeX; ++i) {
    if (i * dx < 0.5) {
      DensP[i] = RHOL;
      Pres[i] = PL;
      Xvel[i] = XVELL;
      Yvel[i] = YVELL;
      Zvel[i] = ZVELL;
    } else {
      DensP[i] = RHOR;
      Pres[i] = PR;
      Xvel[i] = XVELR;
      Yvel[i] = YVELR;
      Zvel[i] = ZVELR;
    }
  }
  Prims2Cons(Prims, Cons, 0, REdgeX);
}

void Domain::ShuOsherIC() {

  for (int i = 0; i < REdgeX; ++i) {
    if (dx * (i - XStart) <= 0.5) {
      Dens[i] = 3.857143;
      Xvel[i] = 2.629369;
#if NDIMS > 1
      Yvel[i] u = 0.0;
#endif
      Pres[i] = 10.33333;
    } else {
      Dens[i] =
          1.0 + (0.2 * std::sin(5.0 * (-4.5 + (dx * 0.5 + dx * (i - XStart)))));
      Xvel[i] = 0.;
#if NDIMS > 1
      Yvel[i] = 0.;
#endif
      Pres[i] = 1.;
    }
  }
  Prims2Cons(Prims, Cons, 0, REdgeX);
}
