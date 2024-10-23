#include "../include/DomainClass.hpp"
#include <cblas.h>

void Domain::ShuOsherIC() {

  for (int i = 0; i < REdgeX; ++i) {
    if (dx * (i - XStart) <= 0.5) {
      DENS[i] = 3.857143;
      XVEL[i] = 2.629369;
      PRES[i] = 10.33333;

    } else {
      DENS[i] =
          1.0 + (0.2 * std::sin(5.0 * (-4.5 + (dx * 0.5 + dx * (i - XStart)))));
      XVEL[i] = 0.;
      PRES[i] = 1.;
    }
  }

  Prims2Cons();
}
void Domain::SodIC() {

  for (int i = 0; i < REdgeX; ++i) {
    if (dx * (i - XStart) <= 0.5) {

      DENS[i] = 1.0;
      XVEL[i] = 0.0;
      PRES[i] = 1.0;
    } else {

      DENS[i] = .125;
      XVEL[i] = 0.0;
      PRES[i] = .1;
    }
  }

  Prims2Cons();
}
