#include "../include/DomainClass.hpp"
#include "../include/GenVarConvert.hpp"
#include "../include/Parameters.h"
// #include <cfenv>
#include <iostream>

int main() {
  // feenableexcept(FE_INVALID);

  Domain Solver;

  (Solver.*(Solver.IC))();

  // for (int i = 0; i < REdgeX * NumVar; i++) {
  //   std::cout << Solver.Prims[i] << std::endl;
  // }
  // std::cout << "\n \n" << std::endl;
  // Solver.Prims2Cons();
  // Solver.Cons2Prim();
  // for (int i = 0; i < REdgeX * NumVar; i++) {
  //   std::cout << Solver.Prims[i] << std::endl;
  // }
  // exit(0);

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
