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

  // for (int i = Solver.XStart; i < Solver.REdgeX; ++i) {
  //   std::cout << Solver.DENS[i * Solver.yDim] << " ";
  // }
  (Solver.*(Solver.BC))("Cons");

  std::cout << "BC DIDNT CRASH!" << std::endl;

  return 0;
}
