#include "../include/Flux.hpp"

void Domain::Calculate_Quad_Points() {
  // std::cout << REdgeX << std::endl;

  for (int i = 2; i < REdgeX - 2; ++i) {
    for (int j = 2; j < REdgeY - 2; ++j) {
      (*GetCell(i, j)).GP_Reconstruction();
      // std::cout << GetCell(i, j)->x << " " << i << " " << " "
      //           << GetCell(i, j)->y << std::endl;
    }
  }
  // exit(0);
}

void Domain::SolveRiemann() {
  for (int i = XStart - 1; i < XEnd + 1; ++i) {
    for (int j = YStart - 1; j < YEnd + 1; ++j) {
      (*GetCell(i, j)).RiemannSolver();
    }
  }
}

void Domain::TimeStep() {
  for (int i = XStart; i < XEnd; ++i) {
    for (int j = YStart; j < YEnd; ++j) {
      (*GetCell(i, j)).FluxRecon();
      (*GetCell(i, j)).AcceptUpdate();
    }
  }
}
