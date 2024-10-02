#include "../include/DomainClass.hpp"
#include <cblas.h>

void Domain::ShockTubeIC() {
  for (int i = 0; i < REdgeX; ++i) {
    if ((float)i / (float)REdgeX <= 0.5) {
      DENSP[i] = RhoL;
      PRES[i] = PL;
      XVEL[i] = XVelL;
      YVEL[i] = YVelL;
      ZVEL[i] = ZVelL;
    } else {
      DENSP[i] = RhoR;
      PRES[i] = PR;
      XVEL[i] = XVelR;
      YVEL[i] = YVelR;
      ZVEL[i] = ZVelR;
    }
  }
  Prims2Cons();
}

void Domain::ShuOsherIC() {

  for (int i = 0; i < REdgeX; ++i) {
    if (dx * (i - XStart) <= 0.5) {
      DENS[i] = 3.857143;
      XVEL[i] = 2.629369;
#if NDIMS > 1
      YVEL[i] u = 0.0;
#endif
      PRES[i] = 10.33333;
    } else {
      DENS[i] =
          1.0 + (0.2 * std::sin(5.0 * (-4.5 + (dx * 0.5 + dx * (i - XStart)))));
      XVEL[i] = 0.;
#if NDIMS > 1
      YVEL[i] = 0.;
#endif
      PRES[i] = 1.;
    }
  }
  Prims2Cons();
}
