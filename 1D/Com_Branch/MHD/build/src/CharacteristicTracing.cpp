#include "../include/DomainClass.hpp"

void Domain::CharTracing() {

  Cons2Prim(Cons, Prims, 0, REdgeX);

  (*this.*SpaceRecon)(XStart - 1, XEnd + 2);

  (*this.*RiemannSolver)(XStart - 1, XEnd + 1);

  Recon(XStart, XEnd);

  (*this.*BC)();
}
