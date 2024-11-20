#include "../include/DomainClass.hpp"
#include "../include/Parameters.h"
#include <cfenv>
#include <iostream>

int main() {
  feenableexcept(FE_INVALID);

  Domain Solver;

  (Solver.*(Solver.IC))();
  (Solver.*(Solver.BC))();

  // Solver.Prims[Tidx(0, 10)] = 10.0;
  // Solver.Prims[Tidx(1, 10)] = .3;
  // Solver.Prims[Tidx(2, 10)] = .4;
  // Solver.Prims[Tidx(3, 10)] = .1;
  // Solver.Prims[Tidx(4, 10)] = 30.0;

  // for (int var = 0; var < NumVar; ++var) {
  //   std::cout << Solver.Prims[Tidx(var, 10)] << " ";
  // }
  // std::cout << std::endl;

  // std::cout << std::endl;
  // Solver.Find_Cs(Solver.Prims, Solver.Cs, 10, 11);
  // Solver.Chars.EigenVectors(Solver.Prims, Solver.Cs, 10, 11);
  // Solver.Chars.Prim2Char(Solver.Prims, 10, 11);

  // Solver.Chars.DumpEigVec(10);
  // std::cout << std::endl;
  // Solver.Chars.DumpChars(10, 11);

  // std::cout << std::endl;
  // Solver.Chars.Char2Prim(Solver.Chars.w, Solver.Prims, 10, 11);

  // for (int var = 0; var < NumVar; ++var) {
  //   std::cout << Solver.Prims[Tidx(var, 10)] << " ";
  // }

  // exit(0);
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
