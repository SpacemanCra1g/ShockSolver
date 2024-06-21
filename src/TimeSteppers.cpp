#include "../include/DomainClass.hpp"

void Domain::ForwardEuler() {

  Calculate_Quad_Points();

  Flux_Recon();
}

void Domain::RK3() {

  for (int i = 0; i < 4; ++i) {
    std::copy(Cons[i], Cons[i] + xDim * yDim, CopyCons[i]);
  }

  ForwardEuler;

  BC

      ForwardEuler;

  Domain_ADD;

  Forward Euler;

  Domain_Add
}
