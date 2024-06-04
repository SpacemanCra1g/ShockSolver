#include "../include/DomainClass.hpp"
#include <cblas.h>
#include <cmath>

void Domain::NeumannBC(std::string Vars) {
  if (Vars == "Cons") {
    for (int i = 0; i < ngc; ++i) {
      // Copies the Left edge of the X dimention
      cblas_dcopy(yDim, &DENS[XStart * yDim], 1, &DENS[yDim * i], 1);
      cblas_dcopy(yDim, &ENERGY[XStart * yDim], 1, &ENERGY[yDim * i], 1);
      cblas_dcopy(yDim, &MOMX[XStart * yDim], 1, &MOMX[yDim * i], 1);
      cblas_dcopy(yDim, &MOMY[XStart * yDim], 1, &MOMY[yDim * i], 1);

      // Copies the Right edge of the X dimention
      cblas_dcopy(yDim, &DENS[(XEnd - 1) * yDim], 1, &DENS[(XEnd + i) * yDim],
                  1);
      cblas_dcopy(yDim, &ENERGY[(XEnd - 1) * yDim], 1,
                  &ENERGY[(XEnd + i) * yDim], 1);
      cblas_dcopy(yDim, &MOMX[(XEnd - 1) * yDim], 1, &MOMX[(XEnd + i) * yDim],
                  1);
      cblas_dcopy(yDim, &MOMY[(XEnd - 1) * yDim], 1, &MOMY[(XEnd + i) * yDim],
                  1);

      // Copies the Left edge of the Y dimention
      cblas_dcopy(xDim, &DENS[YStart], yDim, &DENS[i], yDim);
      cblas_dcopy(xDim, &ENERGY[YStart], yDim, &ENERGY[i], yDim);
      cblas_dcopy(xDim, &MOMX[YStart], yDim, &MOMX[i], yDim);
      cblas_dcopy(xDim, &MOMY[YStart], yDim, &MOMY[i], yDim);

      // Copies the Right edge of the Y dimention
      cblas_dcopy(xDim, &DENS[YEnd - 1], yDim, &DENS[YEnd + i], yDim);
      cblas_dcopy(xDim, &ENERGY[YEnd - 1], yDim, &ENERGY[YEnd + i], yDim);
      cblas_dcopy(xDim, &MOMX[YEnd - 1], yDim, &MOMX[YEnd + i], yDim);
      cblas_dcopy(xDim, &MOMY[YEnd - 1], yDim, &MOMY[YEnd + i], yDim);
    }

  } else if (Vars == "Prims") {
    for (int i = 0; i < ngc; ++i) {
      // Copies the Left edge of the X dimention
      cblas_dcopy(yDim, &DENS[XStart * yDim], 1, &DENS[yDim * i], 1);
      cblas_dcopy(yDim, &PRES[XStart * yDim], 1, &PRES[yDim * i], 1);
      cblas_dcopy(yDim, &XVEL[XStart * yDim], 1, &XVEL[yDim * i], 1);
      cblas_dcopy(yDim, &YVEL[XStart * yDim], 1, &YVEL[yDim * i], 1);

      // Copies the Right edge of the X dimention
      cblas_dcopy(yDim, &DENS[(XEnd - 1) * yDim], 1, &DENS[(XEnd + i) * yDim],
                  1);
      cblas_dcopy(yDim, &PRES[(XEnd - 1) * yDim], 1, &PRES[(XEnd + i) * yDim],
                  1);
      cblas_dcopy(yDim, &XVEL[(XEnd - 1) * yDim], 1, &XVEL[(XEnd + i) * yDim],
                  1);
      cblas_dcopy(yDim, &YVEL[(XEnd - 1) * yDim], 1, &YVEL[(XEnd + i) * yDim],
                  1);

      // Copies the Left edge of the Y dimention
      cblas_dcopy(xDim, &DENS[YStart], yDim, &DENS[i], yDim);
      cblas_dcopy(xDim, &PRES[YStart], yDim, &PRES[i], yDim);
      cblas_dcopy(xDim, &XVEL[YStart], yDim, &XVEL[i], yDim);
      cblas_dcopy(xDim, &YVEL[YStart], yDim, &YVEL[i], yDim);

      // Copies the Right edge of the Y dimention
      cblas_dcopy(xDim, &DENS[YEnd - 1], yDim, &DENS[YEnd + i], yDim);
      cblas_dcopy(xDim, &PRES[YEnd - 1], yDim, &PRES[YEnd + i], yDim);
      cblas_dcopy(xDim, &XVEL[YEnd - 1], yDim, &XVEL[YEnd + i], yDim);
      cblas_dcopy(xDim, &YVEL[YEnd - 1], yDim, &YVEL[YEnd + i], yDim);
    }
  }

  else {
    std::cout << "Invalid call to ApplyBC, terminating" << std::endl;
    exit(1);
  }
}
void Domain::ShuOsherBC(std::string Vars) {
  if (Vars == "Cons") {
    for (int i = 0; i < ngc; ++i) {

      // Copies the Far Left edge of the X dimention
      // This should never change from the inits, so should be good.
      cblas_dcopy(yDim, &DENS[0], 1, &DENS[yDim * i], 1);
      cblas_dcopy(yDim, &ENERGY[0], 1, &ENERGY[yDim * i], 1);
      cblas_dcopy(yDim, &MOMX[0], 1, &MOMX[yDim * i], 1);
      cblas_dcopy(yDim, &MOMY[0], 1, &MOMY[yDim * i], 1);

      // Copies the Right edge of the X dimention
      // We can do something similar to what we did for the left end
      // We need some slight changes for the DENS variable however

      cblas_dcopy(yDim, &ENERGY[(REdgeX - 1) * yDim], 1,
                  &ENERGY[(XEnd + i) * yDim], 1);

      std::fill(&MOMY[(XEnd + i) * yDim], &MOMY[(XEnd + i + 1) * yDim], 0.0);
      std::fill(&MOMX[(XEnd + i) * yDim], &MOMX[(XEnd + i + 1) * yDim], 0.0);

      cblas_dcopy(yDim, &DENS[(REdgeX - 1) * yDim], 1, &DENS[(XEnd + i) * yDim],
                  1);

      // Copies the Left edge of the Y dimention
      cblas_dcopy(xDim, &DENS[YStart], yDim, &DENS[i], yDim);
      cblas_dcopy(xDim, &ENERGY[YStart], yDim, &ENERGY[i], yDim);
      cblas_dcopy(xDim, &MOMX[YStart], yDim, &MOMX[i], yDim);
      cblas_dcopy(xDim, &MOMY[YStart], yDim, &MOMY[i], yDim);

      // Copies the Right edge of the Y dimention
      cblas_dcopy(xDim, &DENS[YEnd - 1], yDim, &DENS[YEnd + i], yDim);
      cblas_dcopy(xDim, &ENERGY[YEnd - 1], yDim, &ENERGY[YEnd + i], yDim);
      cblas_dcopy(xDim, &MOMX[YEnd - 1], yDim, &MOMX[YEnd + i], yDim);
      cblas_dcopy(xDim, &MOMY[YEnd - 1], yDim, &MOMY[YEnd + i], yDim);
    }
  }
}
