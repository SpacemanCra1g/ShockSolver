#include "../include/FluxClass.hpp"

void FluxClass::Mood(int quad, int var, int x) {
  for (int dir = Left; dir <= Right; dir++) {
    if (std::fmin(MoodOrd[x], MoodOrd[x - 1 + 2 * dir]) == 5) {
      GPR2Side(Cons, quad, var, x, dir);
      // GPR2(quad, var, x, y);
    } else if (std::fmin(MoodOrd[x], MoodOrd[x - 1 + 2 * dir]) == 3) {

      GPR1Side(Cons, quad, var, x, dir);
    } else {

      FOGSide(Cons, quad, var, x, dir);
    }
  }
}
