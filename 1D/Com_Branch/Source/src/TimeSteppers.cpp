#include "../include/DomainClass.hpp"
// #define PStuff

void Domain::ForwardEuler() {
  double *PrintVar;

  Cons2Prim(Cons, Prims, 0, REdgeX);

#if SpaceMethod == MOOD
  std::copy(Prims, Prims + NumVar * xDim, PrimsCopy);
  std::fill(Troubled, Troubled + xDim, true);
#endif

#ifdef PStuff

  PrintVar = Prims;
  std::cout << "Prims" << std::endl;
  for (int i = 0; i < xDim; ++i) {
    for (int var = 0; var < NumVar; ++var) {
      std::cout << PrintVar[Tidx(var, i)] << " ";
    }
    std::cout << "   Cell Number " << i << std::endl;
  }

  std::cout << std::endl;
  std::cout << "Prims copy" << std::endl;
  PrintVar = PrimsCopy;
  for (int i = 0; i < xDim; ++i) {
    for (int var = 0; var < NumVar; ++var) {
      std::cout << PrintVar[Tidx(var, i)] << " ";
    }
    std::cout << "   Cell Number " << i << std::endl;
  }

  std::cout << std::endl;
  std::cout << "Troubled" << std::endl;

  for (int i = 0; i < xDim; ++i) {

    std::cout << Troubled[i] << " ";

    std::cout << "   Cell Number " << i << std::endl;
  }
#endif

  (*this.*SpaceRecon)(XStart - 1, XEnd + 1);

#ifdef PStuff
  std::cout << std::endl;
  std::cout << "Left States" << std::endl;
  PrintVar = FluxWalls_Prims[LEFT];
  for (int i = 0; i < xDim; ++i) {
    for (int var = 0; var < NumVar; ++var) {
      std::cout << PrintVar[Tidx(var, i)] << " ";
    }
    std::cout << "   Cell Number " << i << std::endl;
  }

  std::cout << std::endl;
  std::cout << "Right States" << std::endl;
  PrintVar = FluxWalls_Prims[RIGHT];
  for (int i = 0; i < xDim; ++i) {
    for (int var = 0; var < NumVar; ++var) {
      std::cout << PrintVar[Tidx(var, i)] << " ";
    }
    std::cout << "   Cell Number " << i << std::endl;
  }
  // exit(0);
#endif

  MoodFinished = false;

  (this->*RiemannSolver)(XStart - 1, XEnd);

#ifdef PStuff

  std::cout << std::endl;
  std::cout << "Post Riemann Cell Flux" << std::endl;
  PrintVar = CellFlux;
  for (int i = 0; i < xDim; ++i) {
    for (int var = 0; var < NumVar; ++var) {
      std::cout << PrintVar[Tidx(var, i)] << "    ";
    }
    std::cout << "   Cell Number " << i << std::endl;
  }
  // exit(0);
#endif
  Recon(XStart, XEnd);

#ifdef PStuff
  Cons2Prim(Cons, Prims, 0, xDim);

  std::cout << std::endl;
  std::cout << "Final State after updating" << std::endl;
  PrintVar = Cons;
  for (int i = 0; i < xDim; ++i) {
    for (int var = 0; var < NumVar; ++var) {
      std::cout << PrintVar[Tidx(var, i)] << "         ";
    }
    std::cout << "   Cell Number " << i << std::endl;
  }

#endif

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

#ifdef PStuff
  std::cout << std::endl;
  std::cout << "Troubled" << std::endl;

  for (int i = 0; i < xDim; ++i) {

    std::cout << Troubled[i] << " ";

    std::cout << "   Cell Number " << i << std::endl;
  }
#endif
#endif

  // exit(0);
  (*this.*BC)();
}

void Domain::RK3() {
#if SpaceMethod == MOOD
  std::fill(MoodOrd, MoodOrd + xDim, MoodOrder);
  std::fill(DMP_MaxRho, DMP_MaxRho + xDim, -1.0e14);
  std::fill(DMP_MinRho, DMP_MinRho + xDim, 1.0e14);
  std::fill(U2_MaxC, U2_MaxC + xDim, -1.0e14);
  std::fill(U2_MinC, U2_MinC + xDim, 1.0e14);
#endif

  DomainCopy(Cons, CopyBuffer);

  ForwardEuler();

  ForwardEuler();

  DomainAdd(.75, .25, CopyBuffer, Cons);

  ForwardEuler();
  DomainAdd(1.0 / 3.0, 2.0 / 3.0, CopyBuffer, Cons);
}

// void Domain::RK4() {
//   double c1 = 0.391752226571890;
//   double a20 = 0.444370493651235;
//   double c2 = 0.368410593050371;
//   double a21 = 0.555629506348765;
//   double a30 = 0.620101851488403;
//   double c3 = 0.251891774271694;
//   double a32 = 0.379898148511597;
//   double c4 = 0.544974750228521;
//   double a40 = 0.178079954393132;
//   double a43 = 0.821920045606868;
//   double f4 = 0.386708617503269;
//   double ff4 = 0.226007483236906;
//   double ff3 = 0.063692468666290;
//   double f2 = 0.517231671970585;
//   double f3 = 0.096059710526147;
// #if SpaceMethod == Mood53
//   std::fill(Flux.MoodOrd, Flux.MoodOrd + xDim, 5);
//   std::copy(Cons, Cons + NumVar * xDim, ConsCopy);
// #endif

//   DomainCopy(Cons, CopyBuffer);

//   ForwardEuler();

//   DomainAdd(c1, 1.0, CopyBuffer, Cons);

//   DomainCopy(Cons, U1);

//   ForwardEuler();

//   // DomainAdd(c2, a20, c2, CopyBuffer, Cons);
//   DomainAdd(a21, 1.0, U1, Cons);

//   DomainCopy(Cons, U2);
//   ForwardEuler();

//   DomainAdd(a30, c3, CopyBuffer, Cons);
//   DomainAdd(a32, 1.0, U2, Cons);

//   DomainCopy(Cons, U3);
//   ForwardEuler();
//   DomainCopy(Cons, FU3);

//   DomainAdd(a40, c4, CopyBuffer, Cons);
//   DomainAdd(a43, 1.0, U3, Cons);

//   DomainCopy(Cons, U4);
//   ForwardEuler();

//   DomainAdd(f4, ff4, U4, Cons);
//   DomainAdd(ff3, 1.0, FU3, Cons);
//   DomainAdd(f3, 1.0, U3, Cons);
//   DomainAdd(f2, 1.0, U2, Cons);
// }
