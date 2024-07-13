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

  // std::cout << Solver.GetCell(10, 17)->y << std::endl;
  // exit(0);

  (Solver.*(Solver.IC))();

  (Solver.*(Solver.BC))("Cons");

  Solver.SolvePressure();

  do {
    Solver.Find_dt();

    Solver.T += Solver.dt;

    // Solver.Check();
    Solver.RK3();
    std::cout << "The time is: " << Solver.T << " dt = :" << Solver.dt
              << std::endl;

    // exit(0);
  } while (Solver.T < Solver.TN || false);

  Solver.Cons2Prim();

  Solver.WriteOutSolution();

  return 0;
}
