#include "../include/DomainClass.hpp"
#include "../include/Parameters.h"
#include "../include/SRVarConvert.hpp"
// #include <cfenv>
#include <iostream>

int main() {
  // feenableexcept(FE_INVALID);

  Domain Solver;

  (Solver.*(Solver.IC))();

  (Solver.*(Solver.BC))("Cons");
  int counter = 0;

  do {
    Solver.Find_dt();

    Solver.T += Solver.dt;

    (Solver.*(Solver.RK_TimeStepper))();

    std::cout << "The time is: " << Solver.T;
    std::cout << " dt = :" << Solver.dt << std::endl;
    counter++;

  } while (Solver.T < TN); // x while (counter < 26);

  Solver.writeResults();
  return 0;
}
