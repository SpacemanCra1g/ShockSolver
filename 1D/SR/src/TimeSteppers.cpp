#include "../include/DomainClass.hpp"

void Domain::ForwardEuler() {

  (*this.*SpaceRecon)(XStart - 1, XEnd + 1);

  // for (int var = 0; var < NumVar; ++var) {
  //   for (int x = 0; x < xDim; ++x) {
  //     std::cout << FluxWalls_Cons[LEFT][Tidx(var, x)] << std::endl;
  //     ;
  //   }
  //   std::cout << "\n \n" << std::endl;
  // }
  // exit(0);

  MoodFinished = false;

  Hll(XStart - 1, XEnd);

  Recon(XStart, XEnd);

#if SpaceMethod == Mood53
  while (!MoodFinished) {
    MoodFinished = Flux.Detection(MoodFinished);
  }
#endif

  (*this.*BC)();
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
}
