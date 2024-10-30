#include "../include/DomainClass.hpp"
#include "../include/Parameters.h"
#include <cfenv>
#include <iomanip>
#include <iostream>

int main() {
  feenableexcept(FE_INVALID);
  std::cout << std::setprecision(15);

  Domain Solver;

  (Solver.*(Solver.IC))();

  (Solver.*(Solver.BC))();

  Solver.count = 0;
  do {

    Solver.Find_dt();

    Solver.T += Solver.dt;

    (Solver.*(Solver.RK_TimeStepper))();
    Solver.count++;

    Solver.PrintProgress();
  }

  while (Solver.T < TN);

  Solver.writeResults();
  return 0;

  // Create Merge Conflict: Diogenes Version 1.0
}
