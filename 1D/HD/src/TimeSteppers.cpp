#include "../include/DomainClass.hpp"

void Domain::ForwardEuler() {

  (*this.*SpaceRecon)(XStart - 1, XEnd + 2);

  MoodFinished = false;

  Hll(XStart - 1, XEnd + 1);

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
