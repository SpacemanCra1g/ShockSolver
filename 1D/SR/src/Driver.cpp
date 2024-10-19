#include "../include/DomainClass.hpp"
#include "../include/Parameters.h"
// #include "../include/SRVarConvert.hpp"
#include <cfenv>
#include <iostream>

int main() {
  feenableexcept(FE_INVALID);

  Domain Solver;

  // Solver.FluxWalls_Prims[RIGHT][Tidx(0, 10)] = 100.;
  // Solver.FluxWalls_Prims[RIGHT][Tidx(1, 10)] = .37;
  // Solver.FluxWalls_Prims[RIGHT][Tidx(2, 10)] = .9;
  // Solver.FluxWalls_Prims[RIGHT][Tidx(3, 10)] = .0;
  // Solver.FluxWalls_Prims[RIGHT][Tidx(4, 10)] = .001;

  // Solver.Prims2Cons(Solver.FluxWalls_Prims[RIGHT],
  // Solver.FluxWalls_Cons[RIGHT],
  //                   10, 11);

  // for (int var = 0; var < NumVar; ++var) {
  //   std::cout << Solver.FluxWalls_Cons[RIGHT][Tidx(var, 10)] << std::endl;
  // }

  // std::cout << std::endl;

  // Solver.Cons2Prim(Solver.FluxWalls_Cons[RIGHT],
  // Solver.FluxWalls_Prims[RIGHT],
  //                  10, 11);
  // for (int var = 0; var < NumVar; ++var) {
  //   std::cout << Solver.FluxWalls_Prims[RIGHT][Tidx(var, 10)] << std::endl;
  // }
  // exit(0);

  (Solver.*(Solver.IC))();
  (Solver.*(Solver.BC))();

  // for (int x = 0; x < xDim; ++x) {
  //   std::cout << Solver.Energy[x] << std::endl;
  // }
  // exit(0);p

  do {

    Solver.Find_dt();

    Solver.T += Solver.dt;
    if (Solver.dt < 0.0) {
      std::cout << "dt broke at Time: " << Solver.T << std::endl;
      break;
    }

    (Solver.*(Solver.RK_TimeStepper))();

    std::cout << "The time is: " << Solver.T << " dt = " << Solver.dt
              << std::endl;

  } while (Solver.T < TN);

  Solver.writeResults();
  return 0;
}
