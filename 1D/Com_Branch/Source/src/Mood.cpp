#include "../include/DomainClass.hpp"
#include <iostream>

void Domain::Mood(int start, int stop) {

  int L, R;
  double value;

  for (int i = start; i < stop; ++i) {
    L = (MoodOrd[i] < MoodOrd[i - 1]) ? MoodOrd[i] : MoodOrd[i - 1];
    R = (MoodOrd[i] < MoodOrd[i + 1]) ? MoodOrd[i] : MoodOrd[i + 1];

    switch (L) {

    case 3:
      for (int var = 0; var < NumVar; ++var) {
        value = 0.0;
        for (int j = 0; j < 3; ++j) {
          value += Cons[Tidx(var, i - 1 + j)] * Ker.R1Right[j];
        }
        FluxWalls_Cons[LEFT][Tidx(var, i)] = value;
      }
      break;

    case 5:
      for (int var = 0; var < NumVar; ++var) {
        value = 0.0;
        for (int j = 0; j < 5; ++j) {
          value += Cons[Tidx(var, i - 2 + j)] * Ker.R2Right[j];
        }
        FluxWalls_Cons[LEFT][Tidx(var, i)] = value;
      }
      break;

    default:
      for (int var = 0; var < NumVar; ++var) {
        FluxWalls_Cons[LEFT][Tidx(var, i)] = Cons[Tidx(var, i)];
      }
      break;
    }

    switch (R) {

    case 3:
      for (int var = 0; var < NumVar; ++var) {
        value = 0.0;
        for (int j = 0; j < 3; ++j) {
          value += Cons[Tidx(var, i - 1 + j)] * Ker.R1Left[j];
        }
        FluxWalls_Cons[RIGHT][Tidx(var, i)] = value;
      }
      break;

    case 5:
      for (int var = 0; var < NumVar; ++var) {
        value = 0.0;
        for (int j = 0; j < 5; ++j) {
          value += Cons[Tidx(var, i - 2 + j)] * Ker.R2Left[j];
        }
        FluxWalls_Cons[RIGHT][Tidx(var, i)] = value;
      }
      break;

    default:
      for (int var = 0; var < NumVar; ++var) {
        FluxWalls_Cons[RIGHT][Tidx(var, i)] = Cons[Tidx(var, i)];
      }
      break;
    }
  }

  Cons2Prim(FluxWalls_Cons[LEFT], FluxWalls_Prims[LEFT], start, stop);
  Cons2Prim(FluxWalls_Cons[RIGHT], FluxWalls_Prims[RIGHT], start, stop);
}
