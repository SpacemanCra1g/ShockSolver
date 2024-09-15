#include "../include/FluxClass.hpp"

void FluxClass::Mood(int quad, int var, int x, int y) {
  for (int dir = Left; dir <= Right; dir++) {
    if (std::fmin(MoodOrd[idx(x, y)], MoodOrd[idx(x - 1 + 2 * dir, y)]) == 5) {
      GPR2Side(quad, var, x, y, dir);
    } else if (std::fmin(MoodOrd[idx(x, y)],
                         MoodOrd[idx(x - 1 + 2 * dir, y)]) == 3) {
      GPR1Side(quad, var, x, y, dir);
    } else {
      FOGSide(quad, var, x, y, dir);
    }
  }
}
