#include "../include/AssignMemberFunctionPointers.hpp"
#include "../include/DomainClass.hpp"
#include "../include/OptionsClass.hpp"
#include <iostream>

int main() {

  opt Opts;
  Opts.ReadInits("Parameters");

  Domain Solver(Opts);

  AssignPointers(Solver, Opts);

  Solver.AssignCells();

  (Solver.*(Solver.IC))();

  (Solver.*(Solver.BC))("Cons");

  Solver.SolvePressure();

  do {
    Solver.Find_dt();

    Solver.T += Solver.dt;

    std::cout << "The time is: " << Solver.T << " dt = :" << Solver.dt
              << std::endl;

    // exit(0);
  } while (Solver.T < Solver.TN);

  std::cout << "BC DIDNT CRASH!" << std::endl;

  return 0;
}
