#include "../include/DomainClass.hpp"
#include "../include/SourceParameters.h"
#include <cfenv>
#include <iostream>

int main() {
#ifndef GP_METHOD
  feenableexcept(FE_INVALID);
#endif

  Domain Solver;

  (Solver.*(Solver.IC))();
  (Solver.*(Solver.BC))();

  do {

    Solver.Find_dt();

    Solver.T += Solver.dt;
    if (Solver.dt < 0.0) {
      std::cout << "dt broke at Time: " << Solver.T << std::endl;
      break;
    }

    (Solver.*(Solver.RK_TimeStepper))();

    std::cout << "The time is: " << Solver.T << " dt = " << Solver.dt
              << std::endl;

  } while (Solver.T < TN);

  Solver.writeResults();
  return 0;
}
