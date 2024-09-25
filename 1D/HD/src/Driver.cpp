#include "../include/DomainClass.hpp"
#include "../include/Parameters.h"
// #include <cfenv>
#include <iostream>

int main() {
  // feenableexcept(FE_INVALID);

  Domain Solver;

  Solver.SolutionKer.calculate_Preds1D(2);

  (Solver.*(Solver.IC))();

  (Solver.*(Solver.BC))("Cons");

  do {
    Solver.Find_dt();

    Solver.T += Solver.dt;

    (Solver.*(Solver.RK_TimeStepper))();

    std::cout << "The time is: " << Solver.T;
    std::cout << " dt = :" << Solver.dt << std::endl;

  } while (Solver.T < TN);

  Solver.writeResults();
  return 0;
}
