#include "../include/DomainClass.hpp"

void ReMax(int i, double *P) {

  for (int k = 0; k < i; ++k) {
    if (P[k] < .01) {
      P[k] = .01;
    }
  }
}

void Domain::ForwardEuler() {
  Flux.SpaceRecon();
  MoodFinished = false;

  // for (int var = DensP; var <= Pres; ++var) {
  //   for (int i = 0; i < REdgeX; ++i) {
  //     std::cout << Flux.FluxDir[Left][0][var][i] << " "
  //               << Flux.FluxDir[Right][0][var][i] << std::endl;
  //   }
  //   std::cout << "$$$$$$$$$$$$$$$$$$$$$$$" << "\n " << "the var is" << var
  //             << "\n"
  //             << std::endl;
  // }
  // exit(0);

  Flux.HLL();
  exit(0);

  // for (int var = DensP; var <= Pres; ++var) {
  //   for (int i = 0; i < REdgeX; ++i) {
  //     std::cout << Flux.Flux[0][var][i] << std::endl;
  //   }
  //   std::cout << "$$$$$$$$$$$$$$$$$$$$$$$" << "\n " << "the var is" << var
  //             << "\n"
  //             << std::endl;
  // }
  // exit(0);

  Flux.Recon();

#if SpaceMethod == Mood53
  while (!MoodFinished) {
    MoodFinished = Flux.Detection(MoodFinished);
  }
#endif

  (*this.*BC)("Cons");
}

void Domain::RK3() {
#if SpaceMethod == Mood53
  std::fill(Flux.MoodOrd, Flux.MoodOrd + xDim, 5);
  std::copy(Cons, Cons + NumVar * xDim, ConsCopy);
#endif

  DomainCopy();

  ForwardEuler();

  ForwardEuler();

  DomainAdd(.75, .25);

  ForwardEuler();
  DomainAdd(1.0 / 3.0, 2.0 / 3.0);
  ReMax(REdgeX, &Prims[Tidx(Pres, 0)]);
}
