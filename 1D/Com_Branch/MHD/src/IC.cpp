#include "../include/DomainClass.hpp"
#include <cblas.h>

void Domain::ShockTubeIC() {
  for (int i = 0; i < REdgeX; ++i) {
    if ((i - NGC) * dx + dx * 0.5 < 0.5) {
      DensP[i] = 1.0;
      Pres[i] = 1.0;
      Xvel[i] = 0.0;
      Yvel[i] = 0.0;
      Zvel[i] = 0.0;
      PMagX[i] = 0.0;
      PMagY[i] = 0.0;
      PMagZ[i] = 0.0;
    } else {
      DensP[i] = 0.125;
      Pres[i] = 0.1;
      Xvel[i] = 0.0;
      Yvel[i] = 0.0;
      Zvel[i] = 0.0;
      PMagX[i] = 0.0;
      PMagY[i] = 0.0;
      PMagZ[i] = 0.0;
    }
  }
  Prims2Cons(Prims, Cons, 0, REdgeX);
}

void Domain::BrioWu() {
  for (int i = 0; i < REdgeX; ++i) {
    if ((i - NGC) * dx + dx * 0.5 < 0.5) {
      DensP[i] = 1.0;
      Pres[i] = 1.0;
      Xvel[i] = 0.0;
      Yvel[i] = 0.0;
      Zvel[i] = 0.0;
      PMagX[i] = 0.75;
      PMagY[i] = 1.0;
      PMagZ[i] = 0.0;
    } else {
      DensP[i] = 0.125;
      Pres[i] = 0.1;
      Xvel[i] = 0.0;
      Yvel[i] = 0.0;
      Zvel[i] = 0.0;
      PMagX[i] = 0.75;
      PMagY[i] = -1.0;
      PMagZ[i] = 0.0;
    }
  }
  Prims2Cons(Prims, Cons, 0, REdgeX);
}

void Domain::AlfvenWave() {
  double x;
  for (int i = 0; i < REdgeX; ++i) {
    x = (i - NGC) * dx + dx * 0.5;
    if (x < 0.7 && x > 0.3) {
      DensP[i] = 1.0;
      Pres[i] = 1.0;
      Xvel[i] = 0.0;
      Yvel[i] = .0; //-1.5 * (x - .3) * (x - .7);
      Zvel[i] = 0.0;
      PMagX[i] = 0.8;
      PMagY[i] = -2.5 * (x - 0.3) * (x - .7);
      PMagZ[i] = 0.0;
    } else {
      DensP[i] = 1.0;
      Pres[i] = 1.0;
      Xvel[i] = 0.0;
      Yvel[i] = .0;
      Zvel[i] = 0.0;
      PMagX[i] = 0.8;
      PMagY[i] = 0.0;
      PMagZ[i] = 0.0;
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

void Domain::SlowShockIC() {
  for (int i = 0; i < REdgeX; ++i) {
    if (-15.0 + dx * (i - XStart) <= 0.0) {
      // DensP[i] = 5.99924;
      // Xvel[i] = 19.5975;
      // Pres[i] = 460.894;

      DensP[i] = 5.6698;
      Xvel[i] = -1.4701;
      Pres[i] = 100.0;
    } else {
      // DensP[i] = 5.99242;
      // Xvel[i] = -6.19633;
      // Pres[i] = 46.0950;

      DensP[i] = 1.0;
      Xvel[i] = -10.5;
      Pres[i] = 1.0;
    }
  }
  Prims2Cons(Prims, Cons, 0, REdgeX);
}

void Domain::RarefactionIC() {
  for (int i = 0; i < REdgeX; ++i) {
    if ((i - NGC) * dx + dx * 0.5 < 0.5) {
      DensP[i] = 1.0;
      Pres[i] = .4;
      Xvel[i] = -2.0;
    } else {
      DensP[i] = 1.0;
      Pres[i] = .4;
      Xvel[i] = 2.0;
    }
  }
  Prims2Cons(Prims, Cons, 0, REdgeX);
}
