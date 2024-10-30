#include "../include/DomainClass.hpp"
#include <cblas.h>
#include <cmath>

void Domain::ShuOsherIC() {
  double xSpot;
  for (int i = 0; i < REdgeX; ++i) {
    xSpot = -4.5 + (dx * 0.5 + dx * (i - XStart));
    // if (XDomain[i] - X0 <= 0.5) {
    if (xSpot + 4.5 <= 0.5) {
      DensP[i] = 3.857143;
      Xvel[i] = 2.629369;
      Pres[i] = 10.33333;
    } else {
      DensP[i] = 1.0 + (0.2 * std::sin(5.0 * xSpot));
      Xvel[i] = 0.;
      Pres[i] = 1.;
    }
    XDomain[i] = xSpot;
  }

  Prims2Cons(Prims, Cons, 0, REdgeX);
}
void Domain::ShockTubeIC() {

  for (int i = 0; i < REdgeX; ++i) {
    XDomain[i] = X0 + dx * (i - XStart) + (0.5 * dx);
    if (XDomain[i] <= 0.5) {
      DensP[i] = 1.0;
      Xvel[i] = 0.0;
      Pres[i] = 1.0;
    } else {
      DensP[i] = .125;
      Xvel[i] = 0.0;
      Pres[i] = .1;
    }
  }

  Prims2Cons(Prims, Cons, 0, REdgeX);
}
