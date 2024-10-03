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
  Flux.SpaceRecon();
  MoodFinished = false;

  Flux.HLL();

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
}
