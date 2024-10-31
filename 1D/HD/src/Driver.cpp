#include "../include/DomainClass.hpp"
#include "../include/Parameters.h"
// #include <cfenv>
#include <iomanip>
#include <iostream>

int main() {
  // feenableexcept(FE_INVALID);

  std::cout << std::setprecision(15);

  Domain Solver;

  Solver.SolutionKer.calculate_Preds1D(2);

  (Solver.*(Solver.IC))();

  (Solver.*(Solver.BC))("Cons");

  Solver.count = 0;

  do {
    Solver.Find_dt();

    Solver.T += Solver.dt;

    // std::cout << "dt is " << Solver.dt << std::endl;
    // exit(0);

    (Solver.*(Solver.RK_TimeStepper))();
    Solver.count++;

    // for (int i = 0; i < xDim; ++i) {
    //   std::cout << Solver.XVEL[i] << std::endl;
    // }

    std::cout << "The time is: " << Solver.T;
    std::cout << " dt = :" << Solver.dt << std::endl;

    // if (Solver.count == 1) {

    //   Solver.Cons2Prim();
    //   for (int i = 0; i < xDim; ++i) {
    //     // for (int var = 0; var < NumVar; ++var) {
    //     // std::cout << Solver.Cons[Tidx(var, i)] << " ";
    //     std::cout << Solver.DENS[i] << " " << Solver.XVEL[i] << " "
    //               << Solver.PRES[i];
    //     // }
    //     std::cout << "Cell Number: " << i << std::endl;
    //   }
    //   // exit(0);
    // }

  } while (Solver.T < TN);

  // std::cout << "\n\n STarting dump" << std::endl;
  // for (int i = 0; i < xDim; ++i) {
  //   std::cout << Solver.MOMX[i] << std::endl;
  // }

  Solver.Cons2Prim();

  // std::cout << "\n\n Vel Dump" << std::endl;
  // for (int i = 0; i < xDim; ++i) {
  //   std::cout << Solver.XVEL[i] << std::endl;
  // }

  // std::cout << "\n\n True Vel" << std::endl;
  // for (int i = 0; i < xDim; ++i) {
  //   std::cout << Solver.MOMX[i] / Solver.DENS[i] << std::endl;
  // }

  Solver.writeResults();
  return 0;
}
