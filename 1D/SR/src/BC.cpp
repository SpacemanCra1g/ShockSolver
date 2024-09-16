#include "../include/DomainClass.hpp"
#include <cblas.h>
#include <cmath>

void Domain::NeumannBC(std::string Vars) {
  if (Vars == "Cons") {
    for (int i = 0; i < NGC; ++i) {
      // Copies the Left edge of the X dimention
      cblas_dcopy(1, &DENS[XStart], 1, &DENS[i], 1);
      cblas_dcopy(1, &ENERGY[XStart], 1, &ENERGY[i], 1);
      cblas_dcopy(1, &MOMX[XStart], 1, &MOMX[i], 1);
      cblas_dcopy(1, &MOMY[XStart], 1, &MOMY[i], 1);

      // Copies the Right edge of the X dimention
      cblas_dcopy(1, &DENS[(XEnd - 1)], 1, &DENS[(XEnd + i)], 1);
      cblas_dcopy(1, &ENERGY[(XEnd - 1)], 1, &ENERGY[(XEnd + i)], 1);
      cblas_dcopy(1, &MOMX[(XEnd - 1)], 1, &MOMX[(XEnd + i)], 1);
      cblas_dcopy(1, &MOMY[(XEnd - 1)], 1, &MOMY[(XEnd + i)], 1);
    }

  } else if (Vars == "Prims") {
    for (int i = 0; i < NGC; ++i) {
      // Copies the Left edge of the X dimention
      cblas_dcopy(1, &DENS[XStart], 1, &DENS[i], 1);
      cblas_dcopy(1, &PRES[XStart], 1, &PRES[i], 1);
      cblas_dcopy(1, &XVEL[XStart], 1, &XVEL[i], 1);
      cblas_dcopy(1, &YVEL[XStart], 1, &YVEL[i], 1);

      // Copies the Right edge of the X dimention
      cblas_dcopy(1, &DENS[(XEnd - 1)], 1, &DENS[(XEnd + i)], 1);
      cblas_dcopy(1, &PRES[(XEnd - 1)], 1, &PRES[(XEnd + i)], 1);
      cblas_dcopy(1, &XVEL[(XEnd - 1)], 1, &XVEL[(XEnd + i)], 1);
      cblas_dcopy(1, &YVEL[(XEnd - 1)], 1, &YVEL[(XEnd + i)], 1);
    }
  }

  else {
    std::cout << "Invalid call to ApplyBC, terminating" << std::endl;
    exit(1);
  }
}
void Domain::ShuOsherBC(std::string Vars) {
  if (Vars == "Cons") {
    for (int i = 0; i < NGC; ++i) {

      for (int var = 0; var < NDIMS; ++var) {
        // Left edge of X direction
        std::fill(&Cons[Tidx(var, 1)], &Cons[Tidx(var, NGC)],
                  Cons[Tidx(var, 0)]);

        // Right edge of X direction
        std::fill(&Cons[Tidx(var, XEnd)], &Cons[Tidx(var, REdgeX)],
                  Cons[Tidx(var, REdgeX - 1)]);
      }
    }
  }
}
