#ifndef CHARACTERISTICCLASS_H_
#define CHARACTERISTICCLASS_H_

#include "Parameters.h"
#include "definitions.hpp"
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

  void EigenVectors(double *Prims, double *Cs, int Start, int Stop) {
    int i = 0;
    double vx, vy, vz, cs2, cs, v2tan, eta, d, h, lor, sroot;
    double **LL, **RR, *lam;

    for (i = Start; i < Stop; ++i) {
      vx = Prims[Tidx(VELX, i)];
      vy = Prims[Tidx(VELY, i)];
      vz = Prims[Tidx(VELZ, i)];
      d = Prims[Tidx(DENSP, i)];
      cs2 = Cs[i];
      LL = L[i];
      RR = R[i];
      lam = lambda[i];
      lor = 1.0 / std::sqrt(1.0 - (vx * vx + vy * vy + vz * vz));

#if EOS == IDEAL
      h = 1.0 + (GAMMA / (GAMMA - 1.0)) * Prims[Tidx(PRES, i)] / d;
#endif

      cs = std::sqrt(cs2);
      v2tan = vy * vy + vz * vz;
      eta = std::sqrt(1.0 - vx * vx - cs2 * v2tan);

      for (int ii = 0; ii < 5; ++ii) {
        for (int jj = 0; jj < 5; ++jj) {
          LL[ii][jj] = 0.0;
          RR[ii][jj] = 0.0;
        }
      }

      LL[0][1] = d * lor / (2.0 * cs * eta);
      LL[0][4] = 1.0 / (2.0 * h * cs2);

      LL[1][1] = -d * lor / (2.0 * cs * eta);
      LL[1][4] = 1.0 / (2.0 * h * cs2);

      LL[2][0] = 1.0;
      LL[2][4] = -1.0 / (h * cs2);

      LL[3][1] = vx * vy / (1.0 - vx * vx);
      LL[3][2] = 1.0;
      LL[3][4] = vy / (lor * lor * d * h * (1.0 - vx * vx));

      LL[4][1] = vx * vz / (1.0 - vx * vx);
      LL[4][3] = 1.0;
      LL[4][4] = vz / (lor * lor * d * h * (1.0 - vx * vx));

      RR[0][0] = 1.0;
      RR[0][1] = cs * eta / (d * lor);
      RR[0][2] =
          -cs * vy * (lor * eta * vx + cs) / (lor * lor * d * (1.0 - vx * vx));

      RR[0][3] =
          -cs * vz * (lor * eta * vx + cs) / (lor * lor * d * (1.0 - vx * vx));

      RR[0][4] = h * cs2;

      RR[1][0] = 1.0;
      RR[1][1] = -cs * eta / (d * lor);
      RR[1][2] =
          cs * vy * (lor * eta * vx - cs) / (lor * lor * d * (1.0 - vx * vx));

      RR[1][3] =
          cs * vz * (lor * eta * vx - cs) / (lor * lor * d * (1.0 - vx * vx));

      RR[1][4] = h * cs2;

      RR[2][0] = 1.0;
      RR[3][2] = 1.0;
      RR[4][3] = 1.0;

      sroot = cs2 / (lor * lor * (1 - cs2));
      lam[0] = (vx + std::sqrt(sroot * (1.0 - vx * vx + sroot))) / (1 + sroot);
      lam[1] = (vx - std::sqrt(sroot * (1.0 - vx * vx + sroot))) / (1 + sroot);
      lam[2] = vx;
      lam[3] = vx;
      lam[4] = vx;
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
};

#endif // CHARACTERISTICCLASS_H_
