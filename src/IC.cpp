#include "../include/IC.hpp"
#include "../include/DomainClass.hpp"
#include <cblas.h>

void Domain::ShuOsherIC() {
  for (int i = 0; i < REdgeX; ++i) {
    if (dx / 2.0 + dx * (i - XStart) <= 0.5) {
      DENS[i * yDim] = 3.857143;
      XVEL[i * yDim] = 2.629369;
      YVEL[i * yDim] = 0.0;
      PRES[i * yDim] = 10.33333;
    } else {
      // std::cout << "got here \n" << std::endl;
      DENS[i * yDim] =
          1 + 0.2 * std::sin(5.0 * (-4.5 + (dx / 2 + dx * (i - XStart))));
      XVEL[i * yDim] = 0.;
      YVEL[i * yDim] = 0.;
      PRES[i * yDim] = 1.;
    }
  }

  // We copy the x-axis we just wrote to all rows in y
  for (int i = 1; i < REdgeY; ++i) {
    cblas_dcopy(xDim, &DENS[0], yDim, &DENS[i], yDim);
    cblas_dcopy(xDim, &PRES[0], yDim, &PRES[i], yDim);
    cblas_dcopy(xDim, &XVEL[0], yDim, &XVEL[i], yDim);
    cblas_dcopy(xDim, &YVEL[0], yDim, &YVEL[i], yDim);
  }
  Prims2Cons();
}
