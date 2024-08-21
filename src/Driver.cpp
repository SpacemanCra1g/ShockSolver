#include "../include/AssignMemberFunctionPointers.hpp"
#include "../include/DomainClass.hpp"
#include "../include/Parameters.h"
#include <iostream>

int main() {

  // opt Opts;
  // Opts.ReadInits("Parameters");

  Domain Solver;

  Solver.AssignCells();

  (Solver.*(Solver.IC))();

  (Solver.*(Solver.BC))("Cons");

  Solver.SolvePressure();

  do {
    Solver.Find_dt();

    Solver.T += Solver.dt;

    (Solver.*(Solver.RK_TimeStepper))();

    std::cout << "The time is: " << Solver.T << " dt = :" << Solver.dt
              << std::endl;

    // exit(0);
  } while (Solver.T < TN);

  std::cout << "BC DIDNT CRASH!" << std::endl;

  return 0;
}
