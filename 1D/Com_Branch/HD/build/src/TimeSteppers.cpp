#include "../include/DomainClass.hpp"
// #define PStuff

void Domain::ForwardEuler() {
  // double *PrintVar;

  Cons2Prim(Cons, Prims, 0, REdgeX);

  // for (int i = 0; i < REdgeX; ++i) {
  //   if (Prims[Tidx(DENS, i)] < 0.0) {
  //     Prims[Tidx(DENS, i)] = .01;
  //   }

  //   if (Prims[Tidx(PRES, i)] < 0.0) {
  //     Prims[Tidx(PRES, i)] = .01;
  //   }
  // }

#if SpaceMethod == MOOD
  std::copy(Prims, Prims + NumVar * xDim, PrimsCopy);
  std::fill(Troubled, Troubled + xDim, true);
#endif

  // #ifdef PStuff

  //   PrintVar = Prims;
  //   std::cout << "Prims" << std::endl;
  //   for (int i = 0; i < xDim; ++i) {
  //     for (int var = 0; var < NumVar; ++var) {
  //       std::cout << PrintVar[Tidx(var, i)] << " ";
  //     }
  //     std::cout << "   Cell Number " << i << std::endl;
  //   }

  //   std::cout << std::endl;
  //   std::cout << "Prims copy" << std::endl;
  //   PrintVar = PrimsCopy;
  //   for (int i = 0; i < xDim; ++i) {
  //     for (int var = 0; var < NumVar; ++var) {
  //       std::cout << PrintVar[Tidx(var, i)] << " ";
  //     }
  //     std::cout << "   Cell Number " << i << std::endl;
  //   }

  //   std::cout << std::endl;
  //   std::cout << "Troubled" << std::endl;

  //   for (int i = 0; i < xDim; ++i) {

  //     std::cout << Troubled[i] << " ";

  //     std::cout << "   Cell Number " << i << std::endl;
  //   }
  // #endif

  (*this.*SpaceRecon)(XStart - 1, XEnd + 2);

  // #ifdef PStuff
  //   std::cout << std::endl;
  //   std::cout << "Left States" << std::endl;
  //   PrintVar = FluxWalls_Prims[LEFT];
  //   for (int i = 0; i < xDim; ++i) {
  //     for (int var = 0; var < NumVar; ++var) {
  //       std::cout << PrintVar[Tidx(var, i)] << " ";
  //     }
  //     std::cout << "   Cell Number " << i << std::endl;
  //   }

  //   std::cout << std::endl;
  //   std::cout << "Right States" << std::endl;
  //   PrintVar = FluxWalls_Prims[RIGHT];
  //   for (int i = 0; i < xDim; ++i) {
  //     for (int var = 0; var < NumVar; ++var) {
  //       std::cout << PrintVar[Tidx(var, i)] << " ";
  //     }
  //     std::cout << "   Cell Number " << i << std::endl;
  //   }
  //   // exit(0);
  // #endif

  MoodFinished = false;

  (this->*RiemannSolver)(XStart-1, XEnd);

  // #ifdef PStuff

  //   std::cout << std::endl;
  //   std::cout << "Post Riemann Cell Flux" << std::endl;
  //   PrintVar = CellFlux;
  //   for (int i = 0; i < xDim; ++i) {
  //     for (int var = 0; var < NumVar; ++var) {
  //       std::cout << PrintVar[Tidx(var, i)] << "    ";
  //     }
  //     std::cout << "   Cell Number " << i << std::endl;
  //   }
  //   // exit(0);
  // #endif
  #if RIEMANN != RCM
  Recon(XStart, XEnd);
  #endif

  // #ifdef PStuff
  //   Cons2Prim(Cons, Prims, 0, xDim);

  //   std::cout << std::endl;
  //   std::cout << "Final State after updating" << std::endl;
  //   PrintVar = Cons;
  //   for (int i = 0; i < xDim; ++i) {
  //     for (int var = 0; var < NumVar; ++var) {
  //       std::cout << PrintVar[Tidx(var, i)] << "         ";
  //     }
  //     std::cout << "   Cell Number " << i << std::endl;
  //   }

  // #endif

  // for (int i = 0; i < REdgeX; ++i) {
  //   for (int var = 0; var < NumVar; ++var) {
  //     std::cout << Cons[Tidx(var, i)] << " ";
  //   }
  //   std::cout << std::endl;
  // }
  // exit(0);

#if SpaceMethod == MOOD
  while (!MoodFinished) {
    IdxStop = 0;
    MoodFinished = Detection();
  }

// #ifdef PStuff
//   std::cout << std::endl;
//   std::cout << "Troubled" << std::endl;

//   for (int i = 0; i < xDim; ++i) {

//     std::cout << Troubled[i] << " ";

//     std::cout << "   Cell Number " << i << std::endl;
//   }
// #endif
#endif

  // exit(0);
  (*this.*BC)();
}
void Domain::RK2() {
  DomainCopy(Cons, CopyBuffer);

  ForwardEuler();
  ForwardEuler();

  DomainAdd(.5, .5, CopyBuffer, Cons);
}

void Domain::RK3() {

  DomainCopy(Cons, CopyBuffer);

  ForwardEuler();

  ForwardEuler();

  DomainAdd(.75, .25, CopyBuffer, Cons);

  ForwardEuler();
  DomainAdd(1.0 / 3.0, 2.0 / 3.0, CopyBuffer, Cons);
}

void Domain::RK4() {
  double a10 = 1., c1 = 0.391752226571890;
  double a20 = 0.444370493651235, a21 = 0.555629506348765, c2 = 0.368410593050371;
  double a30 = 0.620101851488403, a32 = 0.379898148511597, c3 = 0.251891774271694;
  double a40 = 0.178079954393132, a43 = 0.821920045606868, c4 = 0.544974750228521;
  double f2  = 0.517231671970585, f3  = 0.096059710526147, f4 = 0.386708617503269;
  double ff3 = 0.063692468666290, ff4 = 0.226007483236906;
  int iter = NumVar*xDim;

  // DomainCopy(Cons, CopyBuffer);
  // ForwardEuler(); // Cons = Fl, CopyBuffer = U0 
  // DomainAdd(a10, c1, CopyBuffer, Cons);
  // DomainCopy(Cons, U1);                                                                                               // U1 and CopyBuffer, and Cons Fine
  // ForwardEuler(); // Cons = Fl, CopyBuffer = U0, U1 = U1
  // DomainAdd(a21, c2, U1, Cons);
  // DomainAdd(a20, 1.0, CopyBuffer, Cons); // Cons = U2, CopyBuffer = U0, U1 = U1
  // DomainCopy(Cons, U1);  // Cons = U2, CopyBuffer = U0, U1 = U2
  // ForwardEuler();  // Cons = Fl, CopyBuffer = U0, U1 = U2
  // DomainCopy(U1, UNew); // Cons = Fl, CopyBuffer = U0, U1 = U2, UNew = U2                                            // U1, CopyBuffer, Cons, UNew
  // DomainAdd(c3, a32, Cons, U1);
  // DomainAdd(a30, 1.0, CopyBuffer, U1); // Cons = Fl, CopyBuffer = U0, U1 = U3, UNew = U2
  // DomainCopy(U1, Cons); // Cons = U3, CopyBuffer = U0, U1 = U3, UNew = U2
  // DomainAdd(f3, f2, U1 , UNew); // Cons = U3, CopyBuffer = U0, U1 = U3, UNew = f2U2 + f3U3
  // ForwardEuler();  // Cons = FU3, CopyBuffer = U0, U1 = U3, UNew = f2U2 + f3U3
  // DomainAdd(ff3, 1.0, Cons, UNew); // Cons = FU3, CopyBuffer = U0, U1 = U3, UNew = f2*U2 + f3*U3 + ff3*FU3
  // DomainAdd(a43, c4, U1, Cons);
  // DomainAdd(a40, 1.0, CopyBuffer, Cons); // Cons = U4, CopyBuffer = U0, U1 = U3, UNew = f2*U2 + f3*U3 + ff3*FU3
  // DomainCopy(Cons, U1); // Cons = U4, CopyBuffer = U0, U1 = U4, UNew = f2*U2 + f3*U3 + ff3*FU3
  // DomainAdd(f4, 1.0, Cons, UNew); // Cons = U4, CopyBuffer = U0, U1 = U4, UNew = f2*U2 + f3*U3 + ff3*FU3 + f4*U4
  // ForwardEuler(); // Cons = FU4, CopyBuffer = U0, U1 = U4, UNew = f2*U2 + f3*U3 + ff3*FU3 + f4*U4
  // DomainAdd(1.0, ff4, UNew, Cons);

  // std::cout << "Made it to the jump" << std::endl;
  // // cblas_dscal(NumVar * xDim,   0.39681668417970806, Cons, 1);
  // cblas_dscal(NumVar * xDim,   f4, Cons, 1);
  // // cblas_dscal(NumVar * xDim,   .5, Cons, 1);
  // std::cout << "Post Jump" << std::endl;

  // DomainCopy(Cons, Uin);
  for (int i = 0; i < iter; ++i){
    Uin[i] = Cons[i];
  }


  ForwardEuler();
  for (int i = 0; i < iter; ++i){
    U1[i] = Uin[i] + c1*Cons[i];
  }
  // DomainCopy(U1, Cons);
  for (int i = 0; i < iter; ++i){
    Cons[i] = U1[i];
  }

  ForwardEuler();
  for (int i = 0; i < iter; ++i){
    U2[i] = a20*Uin[i] + a21*U1[i] + c2*Cons[i];
  }
  // DomainCopy(U2, Cons);
  for (int i = 0; i < iter; ++i){
    Cons[i] = U2[i];
  }


  ForwardEuler();
  for (int i = 0; i < iter; ++i){
    U3[i] = a30*Uin[i] + a32*U2[i] + c3*Cons[i];
  }
  // DomainCopy(U3, Cons);
  for (int i = 0; i < iter; ++i){
    Cons[i] = U3[i];
  }

  ForwardEuler();
  // DomainCopy(Cons, FU3);
  for (int i = 0; i < iter; ++i){
    FU3[i] = Cons[i];
  }

  for (int i = 0; i < iter; ++i){
    U4[i] = a40*Uin[i] + a43*U3[i] + c4*FU3[i];
  }
  // DomainCopy(U4, Cons);
  for (int i = 0; i < iter; ++i){
    Cons[i] = U4[i];
  }
  ForwardEuler();
  // DomainCopy(Cons, FU4);
  for (int i = 0; i < iter; ++i){
    FU4[i] = Cons[i];
  }
  for (int i = 0; i < iter; ++i){
    Cons[i] = f2*U2[i] + f3*U3[i] + ff3*FU3[i] + f4*U4[i]  + ff4*FU4[i];
  }
}
