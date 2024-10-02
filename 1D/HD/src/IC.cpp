#include "../include/DomainClass.hpp"
#include <cblas.h>

void Domain::ShuOsherIC() {

  for (int i = 0; i < REdgeX; ++i) {
    if (dx * (i - XStart) <= 0.5) {
      //       DENS[i] = 3.857143;
      //       XVEL[i] = 2.629369;
      // #if NDIMS > 1
      //       YVEL[i] u = 0.0;
      // #endif
      //       PRES[i] = 10.33333;
      //
      DENS[i] = 1.0;
      XVEL[i] = 0.0;
      PRES[i] = 1.0;
    } else {
      //       DENS[i] =
      //           1.0 + (0.2 * std::sin(5.0 * (-4.5 + (dx * 0.5 + dx * (i -
      //           XStart)))));
      //       XVEL[i] = 0.;
      // #if NDIMS > 1
      //       YVEL[i] = 0.;
      // #endif
      //       PRES[i] = 1.;
      DENS[i] = .125;
      XVEL[i] = 0.0;
      PRES[i] = .1;
    }
  }

  // We copy the x-axis we just wrote to all rows in y
  // for (int i = 1; i < REdgeY; ++i) {
  //   cblas_dcopy(xDim, &DENS[0], yDim, &DENS[i], yDim);
  //   cblas_dcopy(xDim, &PRES[0], yDim, &PRES[i], yDim);
  //   cblas_dcopy(xDim, &XVEL[0], yDim, &XVEL[i], yDim);
  //   cblas_dcopy(xDim, &YVEL[0], yDim, &YVEL[i], yDim);
  // }
  Prims2Cons();
}
