#include "../include/DomainClass.hpp"
#include "../include/Parameters.h"
#include "../include/SRVarConvert.hpp"
#include <cfenv>
#include <iostream>

int main() {
  feenableexcept(FE_INVALID);

  double maxx = 0.0;
  Domain Solver;

  (Solver.*(Solver.IC))();

  (Solver.*(Solver.BC))("Cons");
  int counter = 0;

  do {
    maxx = 0.0;
    Solver.Find_dt();

    Solver.T += Solver.dt;

    (Solver.*(Solver.RK_TimeStepper))();

    std::cout << "The time is: " << Solver.T << " dt = " << Solver.dt
              << std::endl;

    for (int i = 0; i < REdgeX; ++i) {
      if (maxx < Solver.Prims[Tidx(VelY, i)]) {
        maxx = Solver.Prims[Tidx(VelY, i)];
      }
    }

    counter++;

  } while (Solver.T < TN);
  // } while (counter < 20);

  Solver.writeResults();
  return 0;
}
