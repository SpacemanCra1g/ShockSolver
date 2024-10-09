#include "../include/DomainClass.hpp"
#include "../include/Parameters.h"
#include "../include/SRVarConvert.hpp"
// #include <cfenv>
#include <iostream>

int main() {
  // feenableexcept(FE_INVALID);

  double maxx = 0.0;
  Domain Solver;

  (Solver.*(Solver.IC))();
  // double Lor;
  // double h;
  // double theta;
  // double Test[REdgeX];
  // for (int i = 0; i < REdgeX; ++i) {
  //   for (int var = 0; var < NumVar; ++var) {
  //     // Cons[Tidx(var, i)] = C[var];
  // if (i * dx < 0.5) {
  //   theta = PL / RhoL;
  //   Lor =
  //       std::pow(1.0 - (XVelL * XVelL + YVelL * YVelL + ZVelL * ZVelL), -.5);

  //   h = 1.0 + (GAMMA * theta) / (GAMMA - 1);
  //   Test[i] = Lor * RhoL * Lor * h - PL;

  //       Solver.Cons[Tidx(DensP, i)] = Lor * RhoL;
  //       Solver.Cons[Tidx(MomX, i)] = RhoL * h * Lor * Lor * XVelL;
  //       Solver.Cons[Tidx(MomY, i)] = RhoL * h * Lor * Lor * YVelL;
  //       Solver.Cons[Tidx(MomZ, i)] = RhoL * h * Lor * Lor * ZVelL;
  //       Solver.Cons[Tidx(Ener, i)] = RhoL * h * Lor * Lor - PL;

  // } else {
  //   theta = PR / RhoR;
  //   Lor =
  //       std::pow(1.0 - (XVelR * XVelR + YVelR * YVelR + ZVelR * ZVelR), -.5);

  //   h = 1.0 + (GAMMA * theta) / (GAMMA - 1);

  //   Test[i] = Lor * RhoR * Lor * h - PR;
  //       Solver.Cons[Tidx(MomX, i)] = RhoR * h * Lor * Lor * XVelR;
  //       Solver.Cons[Tidx(MomY, i)] = RhoR * h * Lor * Lor * YVelR;
  //       Solver.Cons[Tidx(MomZ, i)] = RhoR * h * Lor * Lor * ZVelR;
  //       Solver.Cons[Tidx(Ener, i)] = RhoR * h * Lor * Lor - PR;

  // std::cout << Lor * RhoR << " \n"
  //           << RhoR * h * Lor * Lor * XVelR << " \n"
  //           << RhoR * h * Lor * Lor * YVelR << "\n"
  //           << RhoR * h * Lor * Lor * ZVelR << "\n"
  //           << RhoR * h * Lor * Lor - PR << std::endl;
  // exit(0);
  //   }
  //   }
  // }
  // for (int i = 0; i < REdgeX; ++i) {
  //   std::cout << Solver.ENERGY[i] << " " << Test[i] << std::endl;
  // }
  // exit(0);

  (Solver.*(Solver.BC))("Cons");
  int counter = 0;

  do {
    maxx = 0.0;
    Solver.Find_dt();

    Solver.T += Solver.dt;

    (Solver.*(Solver.RK_TimeStepper))();

    std::cout << "The time is: " << Solver.T;

    for (int i = 0; i < REdgeX; ++i) {
      if (maxx < Solver.Prims[Tidx(VelY, i)]) {
        maxx = Solver.Prims[Tidx(VelY, i)];
      }
    }

    std::cout << " The Largest u = " << maxx << std::endl;

    counter++;

  } while (Solver.T < TN); // x while (counter < 26)
  // } while (counter < 1);

  Solver.writeResults();
  return 0;
}
