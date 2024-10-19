#include "../include/DomainClass.hpp"
#include "../include/Parameters.h"

// Most of the structure of this code is taken from the PLUTO solver
// So lots of credit to those authors

#define ENERGY_SOLVE 1
#define PRESURE_FIX_SOLVE 2

void Domain::Prims2Cons(double *Uin, double *Uout, int start, int stop) {
  double v2, Lor, h, alpha;
  for (int i = start; i < stop; ++i) {

    v2 = std::pow(Uin[Tidx(VELX, i)], 2) + std::pow(Uin[Tidx(VELY, i)], 2) +
         std::pow(Uin[Tidx(VELZ, i)], 2);

    Lor = 1.0 / (std::sqrt(1.0 - v2));

#if EOS == IDEAL
    h = 1.0 +
        (GAMMA / (GAMMA - 1.0)) * Uin[Tidx(PRES, i)] / Uin[Tidx(DENSP, i)];
#endif

    alpha = Uin[Tidx(DENSP, i)] * h * Lor * Lor;

    Uout[Tidx(DENS, i)] = Uin[Tidx(DENSP, i)] * Lor;
    Uout[Tidx(MOMX, i)] = Uin[Tidx(VELX, i)] * alpha;
    Uout[Tidx(MOMY, i)] = Uin[Tidx(VELY, i)] * alpha;
    Uout[Tidx(MOMZ, i)] = Uin[Tidx(VELZ, i)] * alpha;
    Uout[Tidx(ENER, i)] = alpha - Uin[Tidx(PRES, i)];
  }
}

void Domain::Cons2Prim(double *Uin, double *Uout, int start, int stop) {
  int SolMethod;
  int err = 0;

  for (int i = start; i < stop; ++i) {
    SolMethod = ENERGY_SOLVE;

    if (SolMethod == ENERGY_SOLVE) {
      err = EnergyInverter(Uin, Uout, i);
      // err = NaiveNewton(Uin, Uout, i);
      if (err) {
        if (err == 1) {
          std::cout << "The equation does not admit a solution" << std::endl;
          std::cout << "Cell Number: " << i << std::endl;
          // double x = std::sqrt(-2.0);
        } else if (err == 2) {
          std::cout << "Negative Pressure" << std::endl;
        } else if (err == 4) {
          std::cout << "Pressure is NaN" << std::endl;
        } else if (err == 3) {
          std::cout << "Other Problem" << std::endl;
        }

        std::cout << "Failure at Time: " << T << std::endl;
        std::cout << "The Cons were :" << std::endl;
        for (int var = 0; var < NumVar; ++var) {
          std::cout << Uin[Tidx(var, i)] << std::endl;
        }
        writeResults();
        exit(0);
        // SolMethod = PRESURE_FIX_SOLVE;
      }
    }

    if (SolMethod == PRESURE_FIX_SOLVE) {
      // err = PressureFix(Uin, Uout, i);
      if (err) {
        std::cout << "CRASH REPORT" << std::endl;
        std::cout << "Failure in Pressure fix" << std::endl;
      }
    }
  }
}

void Domain::Press(int x) {
  double C[NumVar];
  for (int var = 0; var < NumVar; ++var) {
    C[var] = Cons[Tidx(var, x)];
  }
  // PRES[x] = Pressure(C);
}

void Domain::SolvePressure() {
  for (int i = 0; i < xDim; ++i) {
    Press(i);
  }
}

#undef ENERGY_SOLVE
#undef ENTROPY_SOLVE
#undef PRESURE_FIX_SOLVE
