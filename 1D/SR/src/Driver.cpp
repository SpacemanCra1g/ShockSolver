#include "../include/DomainClass.hpp"
#include "../include/Parameters.h"
// #include "../include/SRVarConvert.hpp"
#include <cfenv>
#include <iostream>

int main() {
  feenableexcept(FE_INVALID);

  Domain Solver;

  (Solver.*(Solver.IC))();
  (Solver.*(Solver.BC))();

  do {

    Solver.Find_dt();

    Solver.T += Solver.dt;

    (Solver.*(Solver.RK_TimeStepper))();

    std::cout << "The time is: " << Solver.T << " dt = " << Solver.dt
              << std::endl;

  } while (Solver.T < TN);

  Solver.writeResults();
  return 0;
}
