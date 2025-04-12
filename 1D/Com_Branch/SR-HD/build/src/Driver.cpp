#include "../include/DomainClass.hpp"
#include "../include/SourceParameters.h"
#include <cfenv>
#include <iomanip>
#include <iostream>

int main() {
  // int n = 1;
#ifndef GP_METHOD
  feenableexcept(FE_INVALID);
#endif
  std::cout << std::setprecision(15);
  Domain Solver;

  // for (int j = 0; j < n; ++j) {
  (Solver.*(Solver.IC))();
  (Solver.*(Solver.BC))();
  Solver.T = 0.0;

  // for (int i = 0; i < REdgeX; ++i) {
  //   if ((i - NGC) * dx + dx * 0.5 < 0.5) {
  //     Solver.Yvel[i] += j * (.99 / (double)(n - 1));
  //   } else {
  //     continue;
  //   }
  // }
  // Solver.Prims2Cons(Solver.Prims, Solver.Cons, 0, REdgeX);

  int counter = 0;
  do {
    counter++;
    Solver.Find_dt();

    Solver.T += Solver.dt;
    if (Solver.dt < 0.0) {
      std::cout << "dt broke at Time: " << Solver.T << std::endl;
      break;
    }

    (Solver.*(Solver.RK_TimeStepper))();

    if (counter % 100 == 0) {
      std::cout << "The time is: " << Solver.T << " dt = " << Solver.dt
                << std::endl;
    }
  } while (true && Solver.T < TN);

  Solver.writeResults();
  // }
  return 0;
}
