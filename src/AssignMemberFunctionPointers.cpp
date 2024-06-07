#include "../include/AssignMemberFunctionPointers.hpp"

void AssignPointers(Domain &Solver, opt &Opts) {
  if (Opts.TestProblem == "ShuOsher") {
    Solver.BC = &Domain::ShuOsherBC;
    Solver.IC = &Domain::ShuOsherIC;
  } else {
    if (Opts.BC == "Neumann") {
      Solver.BC = &Domain::NeumannBC;
    } else {
      std::cout << "Invalid Boundary Conditions \nExiting" << std::endl;
      exit(0);
    }
  }

  if(Opts.RK == 3){
    Solver.RK_TimeStepper = &Domain::RK3;
  }
}
