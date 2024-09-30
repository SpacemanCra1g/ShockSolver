#include "../include/DomainClass.hpp"
#include "../include/Parameters.h"
#include "../include/SRVarConvert.hpp"
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
  int counter = 0;

  do {
    Solver.Find_dt();

    Solver.T += Solver.dt;

    (Solver.*(Solver.RK_TimeStepper))();

    std::cout << "The time is: " << Solver.T;
    std::cout << " dt = :" << Solver.dt << std::endl;
    counter++;

  } while (Solver.T < TN); // x while (counter < 26); // (Solver.T < TN);

  // Solver.Cons2Prim();
  // for (int i = 0; i < xDim * NumVar; ++i) {
  //   std::cout << Solver.Prims[i] << std::endl;
  // }

  Solver.writeResults();
  return 0;
}
