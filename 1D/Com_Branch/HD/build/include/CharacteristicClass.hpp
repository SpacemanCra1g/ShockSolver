#ifndef CHARACTERISTICCLASS_H_
#define CHARACTERISTICCLASS_H_

#include "SourceParameters.h"
#include <iostream>

/* Most of the structure of this file and the calculations are inspired heavily
 * by the PLUTO code, specifically Src/RHD/eigenv.c */

class Characteristics {

public:
  double ***L, ***R, **lambda;
  double *w;

  Characteristics() {
    L = new double **[xDim];
    R = new double **[xDim];
    lambda = new double *[xDim];
    // Lambdas = new double *[xDim];
    w = new double[xDim * NumVar];
    for (int i = 0; i < xDim; ++i) {
      lambda[i] = new double[NumVar];
      L[i] = new double *[NumVar];
      R[i] = new double *[NumVar];
      // Lambdas[i] = new double[NumVar];

      for (int j = 0; j < NumVar; ++j) {
        L[i][j] = new double[NumVar];
        R[i][j] = new double[NumVar];
        for (int var = 0; var < NumVar; ++var) {
          L[i][j][var] = 0.0;
          R[i][j][var] = 0.0;
        }
      }
    }
  }

  void Char2Prim(double *Chars, double *Prims, int Start, int Stop) {
    double Val;
    for (int i = Start; i < Stop; ++i) {
      for (int var = 0; var < NumVar; ++var) {
        Val = 0.0;
        for (int waveNum = 0; waveNum < NumVar; ++waveNum) {
          Val += R[i][waveNum][var] * Chars[Tidx(waveNum, i)];
        }
        Prims[Tidx(var, i)] = Val;
      }
    }
  }

  void Char_Recover_Prims(double *State, double *Flux, int Center) {
    double **RR = R[Center];
    double Val;

    for (int var = 0; var < NumVar; ++var) {
      Val = 0.0;
      for (int waveNum = 0; waveNum < NumVar; ++waveNum) {
        Val += RR[waveNum][var] * State[waveNum];
      }
      Flux[Tidx(var, Center)] = Val;
    }
  }

  void Char_Stencil(double **Stencil, int Center, int Radius) {
    double **LL = L[Center];
    double Prim[NumVar];
    for (int vec = 0; vec < 2 * Radius + 1; ++vec) {
      for (int var = 0; var < NumVar; ++var) {
        Prim[var] = Stencil[var][vec];
      }
      for (int waveNum = 0; waveNum < NumVar; ++waveNum) {
        Stencil[waveNum][vec] = 0.0;
        for (int var = 0; var < NumVar; ++var) {
          Stencil[waveNum][vec] += Prim[var] * LL[waveNum][var];
        }
      }
    }
  }

  void DumpChars(int start, int stop) {
    for (int i = start; i < stop; ++i) {
      for (int j = 0; j < NumVar; ++j) {
        std::cout << w[Tidx(j, i)] << " ";
      }
      std::cout << std::endl;
    }
  }

  void DumpEigVec(int i) {
    for (int j = 0; j < NumVar; ++j) {
      for (int k = 0; k < NumVar; ++k) {
        std::cout << L[i][j][k] << " ";
      }
      std::cout << std::endl;
    }
  }

  // Defined in the physics specific module
  // void EigenVectors(double *Prims, double *Cs, int Start, int Stop);
  void EigenVectors(double *Prims, double *Cs, int Start, int Stop);
};

#endif // CHARACTERISTICCLASS_H_
