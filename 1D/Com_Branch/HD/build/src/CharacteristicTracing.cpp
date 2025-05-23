#include "../include/DomainClass.hpp"

void Domain::CharTracing() {

  Cons2Prim(Cons, Prims, 0, REdgeX);

  DomainCopy(Cons, CopyBuffer);
  
  (*this.*SpaceRecon)(XStart - 1, XEnd + 2);

  (*this.*RiemannSolver)(XStart - 1, XEnd + 1);

  Recon(XStart, XEnd);


  #if RIEMANN==HYBRID
  rcm_Counter++;
  Cons2Prim(CopyBuffer,Prims,XStart,XEnd);
  for (int i = XStart; i < XEnd-1; ++i){
    if (RcmReduction[i]){
      for (int var = 0; var < NumVar; ++var){
        FluxWalls_Prims[RIGHT][Tidx(var,i)] = Prims[Tidx(var,i)];
        FluxWalls_Prims[LEFT][Tidx(var,i+1)] = Prims[Tidx(var,i+1)];
        rcm(i-1,i+1);
      }
    }
  }
  Cons2Prim(Cons,Prims,0,REdgeX);
   (*this.*BC)();
  #endif

  (*this.*BC)();
}
