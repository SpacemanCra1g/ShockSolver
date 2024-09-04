#include "../include/DomainClass.hpp"

void Domain::ForwardEuler() {
#if SpaceMethod == Weno
  Flux.WENO();
#elif SpaceMethod == Fog
  Flux.FOG();
#endif
  // Calculate_Quad_Points();
  Flux.HLL();

  Flux.Recon();

  (*this.*BC)("Cons");
}

void Domain::RK3() {

  DomainCopy();

  ForwardEuler();

  ForwardEuler();

  DomainAdd(.75, .25);

  ForwardEuler();
  DomainAdd(1.0 / 3.0, 2.0 / 3.0);
}
