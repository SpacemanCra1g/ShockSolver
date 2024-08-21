#include "../include/DomainClass.hpp"

void Domain::ForwardEuler() {

  Flux.FOG();
  // Calculate_Quad_Points();
  for (int i = 0; i < REdgeX; ++i) {
    std::cout << Flux.RightFlux[0][3][i] << std::endl;
  }
  exit(0);

  // Flux_Recon();
}

void Domain::RK3() {

  // for (int i = 0; i < 4; ++i) {
  //   std::copy(Cons[i], Cons[i] + xDim * yDim, CopyCons[i]);
  // }

  // ForwardEuler;

  // BC

  //     ForwardEuler;

  // Domain_ADD;

  // Forward Euler;

  // Domain_Add
}
