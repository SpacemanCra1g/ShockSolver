#include "../include/DomainClass.hpp"
#include "../include/Parameters.h"
#include <iostream>

int main() {

  Domain Solver;

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
