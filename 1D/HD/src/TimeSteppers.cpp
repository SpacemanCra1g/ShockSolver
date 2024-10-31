#include "../include/DomainClass.hpp"
#include "../include/FluxClass.hpp"

void FluxClass::SpaceRecon() {
  int quad = 0;

  for (int var = 0; var < NumVar; ++var) {
    for (int x = 1; x < REdgeX - 1; ++x) {

#if SpaceMethod == Weno
      WENO(quad, var, x);
#elif SpaceMethod == Fog
      FOG(quad, var, x);
#elif SpaceMethod == Gp1
      GPR1(quad, var, x);
#elif SpaceMethod == Gp2
      GPR2(quad, var, x);
#elif SpaceMethod == Mood53
      Mood(quad, var, x);
#endif
    }
  }
}

void Domain::ForwardEuler() {

  // std::cout << "Inital Cons" << std::endl;
  // for (int i = 0; i < xDim; ++i) {
  //   for (int var = 0; var < NumVar; ++var) {
  //     std::cout << Cons[Tidx(var, i)] << " ";
  //   }
  //   std::cout << "Cell Idx " << i << std::endl;
  // }
  // std::cout << std::endl;

  Flux.SpaceRecon();
  MoodFinished = false;

  Flux.HLL();

  // std::cout << "Flux" << std::endl;
  // for (int i = 0; i < xDim; ++i) {
  //   for (int var = 0; var < NumVar; ++var) {
  //     std::cout << Flux.Flux[0][var][i] << " ";
  //   }
  //   std::cout << "Cell Idx " << i << std::endl;
  // }
  // std::cout << std::endl;

  // for (int i = 0; i < xDim; ++i) {
  //   std::cout << Cons[Tidx(VelX, i)] << std::endl;
  // }

  // std::cout << "\n\n" << std::endl;

  // for (int i = 0; i < xDim; ++i) {
  //   std::cout << Flux.Flux[0][VelX][i] << std::endl;
  // }
  // std::cout << "\n\n" << std::endl;

  Flux.Recon();

  // for (int i = 0; i < xDim; ++i) {
  //   std::cout << Cons[Tidx(VelX, i)] << std::endl;
  // }

  // exit(0);

#if SpaceMethod == Mood53
  while (!MoodFinished) {
    MoodFinished = Flux.Detection(MoodFinished);
  }
#endif

  // std::cout << "Updated Cons" << std::endl;
  // for (int i = 0; i < xDim; ++i) {
  //   for (int var = 0; var < NumVar; ++var) {
  //     std::cout << Cons[Tidx(var, i)] << " ";
  //   }
  //   std::cout << "Cell Idx " << i << std::endl;
  // }
  // std::cout << "\n#####################\n" << std::endl;

  (*this.*BC)("Cons");
}

void Domain::RK3() {
#if SpaceMethod == Mood53
  std::fill(Flux.MoodOrd, Flux.MoodOrd + xDim, 5);
  std::copy(Cons, Cons + NumVar * xDim, ConsCopy);
#endif

  DomainCopy();

  // if (count == 2) {

  // }

  ForwardEuler();



  // if (count == 2) {
  //   std::cout << "Updated Cons" << std::endl;
  //   for (int i = 0; i < xDim; ++i) {
  //     for (int var = 0; var < NumVar; ++var) {
  //       std::cout << Cons[Tidx(var, i)] << " ";
  //     }
  //     std::cout << "Cell Idx " << i << std::endl;
  //   }
  //   std::cout << "\n#####################\n" << std::endl;
  //   std::cout << dt << std::endl;
  //   // exit(0);
  // }

  ForwardEuler();



  DomainAdd(.75, .25);

  ForwardEuler();
  DomainAdd(1.0 / 3.0, 2.0 / 3.0);
// std::cout << "Updated Cons" << std::endl;
//     for (int i = 0; i < xDim; ++i) {
//       for (int var = 0; var < NumVar; ++var) {
//         std::cout << Cons[Tidx(var, i)] << " ";
//       }
//       std::cout << "Cell Idx " << i << std::endl;
//     }
//     std::cout << "\n#####################\n" << std::endl;
//     exit(0);
}
