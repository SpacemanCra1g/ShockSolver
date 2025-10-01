#include "../include/DomainClass.hpp"
#include <cblas.h>

void Domain::ShockTubeIC() {
  for (int i = 0; i < REdgeX; ++i) {
    if ((i - NGC) * dx + dx * 0.5 < 0.5) {
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

void Domain::ShockOnly() {
  for (int i = 0; i < REdgeX; ++i) {
    if ((i - NGC) * dx + dx * 0.5 < 0.5) {
      DensP[i] = .1785866 * 25.;
      Pres[i] = .9040;
      Xvel[i] = .319371;
      Yvel[i] = .77208971;
      Zvel[i] = 0.0;
    } else {
      DensP[i] = .039998 * 25.;
      Pres[i] = .01;
      Xvel[i] = 0.0;
      Yvel[i] = 0.9;
      Zvel[i] = 0.0;
    }
  }
  Prims2Cons(Prims, Cons, 0, REdgeX);
}

void Domain::ContactOnly() {
  for (int i = 0; i < REdgeX; ++i) {
    if ((i - NGC) * dx + dx * 0.5 < 0.5) {
      DensP[i] = .00059660 * 25.;
      Pres[i] = .9046;
      Xvel[i] = .31937058;
      Yvel[i] = .94721706;
      Zvel[i] = 0.0;
    } else {
      DensP[i] = .17858635 * 25.;
      Pres[i] = .9046;
      Xvel[i] = 0.31937058;
      Yvel[i] = 0.772089;
      Zvel[i] = 0.0;
    }
  }
  Prims2Cons(Prims, Cons, 0, REdgeX);
}

void Domain::WenProblem1() {
  for (int i = 0; i < REdgeX; ++i) {
    if ((i - NGC) * dx + dx * 0.5 < 0.5) {
      DensP[i] = 10.0;
      Pres[i] = 13.3;
      Xvel[i] = 0.0;
      Yvel[i] = 0.0;
      Zvel[i] = 0.0;
    } else {
      DensP[i] = 1.0;
      Pres[i] = 0.00000066;
      Xvel[i] = 0.0;
      Yvel[i] = 0.0;
      Zvel[i] = 0.0;
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
