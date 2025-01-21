#include "../include/DomainClass.hpp"

// Prims, CopyBuffer, &Chars, start, stop);

#define sign(a) ((a >= 0.0) ? 1.0 : -1.0)
using namespace std;
void Domain::minmod(double *P, double *dW, Characteristics *Ch, int start,
                    int stop) {

  double Left, Right;
  double *LeftEig;

  for (int i = start; i < stop; ++i) {
    for (int var = 0; var < NumVar; ++var) {
      LeftEig = Ch->L[i][var];
      Left = 0.0;
      Right = 0.0;
      for (int wave = 0; wave < NumVar; ++wave) {
        Left += LeftEig[wave] * (P[Tidx(wave, i + 1)] - P[Tidx(wave, i)]);
        Right += LeftEig[wave] * (P[Tidx(wave, i)] - P[Tidx(wave, i - 1)]);
      }
      dW[Tidx(var, i)] =
          0.5 * (sign(Left) + sign(Right)) * fmin(fabs(Left), fabs(Right));
    }
  }
}

void Domain::mc(double *P, double *dW, Characteristics *Ch, int start,
                int stop) {

  double Left, Right;
  double *LeftEig;

  for (int i = start; i < stop; ++i) {
    for (int var = 0; var < NumVar; ++var) {
      LeftEig = Ch->L[i][var];
      Left = 0.0;
      Right = 0.0;
      for (int wave = 0; wave < NumVar; ++wave) {
        Left += LeftEig[wave] * (P[Tidx(wave, i + 1)] - P[Tidx(wave, i)]);
        Right += LeftEig[wave] * (P[Tidx(wave, i)] - P[Tidx(wave, i - 1)]);
      }
      dW[Tidx(var, i)] =
          (sign(Left) + sign(Right)) *
          fmin(fmin(fabs(Left), fabs(Right)), .25 * fabs(Left + Right));
      // 0.5 * (sign(Left) + sign(Right)) * fmin(fabs(Left), fabs(Right));
    }
  }
}

void Domain::vanleer(double *P, double *dW, Characteristics *Ch, int start,
                     int stop) {
  double Left, Right;
  double *LeftEig;

  for (int i = start; i < stop; ++i) {
    for (int var = 0; var < NumVar; ++var) {
      LeftEig = Ch->L[i][var];
      Left = 0.0;
      Right = 0.0;
      for (int wave = 0; wave < NumVar; ++wave) {
        Left += LeftEig[wave] * (P[Tidx(wave, i + 1)] - P[Tidx(wave, i)]);
        Right += LeftEig[wave] * (P[Tidx(wave, i)] - P[Tidx(wave, i - 1)]);
      }
      dW[Tidx(var, i)] =
          (Left * Right <= 0.0) ? 0.0 : 2.0 * Left * Right / (Left + Right);
    }
  }
}
