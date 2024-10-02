#ifndef FLUXCLASS_H_
#define FLUXCLASS_H_

#include "../include/GP_Kernel.hpp"
#include "../include/Parameters.h"
#include <cmath>

class FluxClass {

public:
  GP_Kernel *Kern;
  double *Cons;
  double *dt;
  bool *Troubled;
  int *MoodOrd;
  bool MoodFinished = true;
  double *Uin;

  double FluxDir[NDIMS * 2][nqp][NumVar][xDim * yDim];

  double Flux[nqp][NumVar][xDim * yDim];

  void Fluxinit(double *COns, double *DT) {
    Cons = COns;
    dt = DT;
  }

  // Defined in the src/FOG.cpp file
  void FOG(int, int, int, int);

  // Defined in the src/GP-FVM.cpp file
  void GP(int, int, int, int);
  void GPR1(int, int, int, int);
  void GPR2(int, int, int, int);
  void GPR1Side(double *, int, int, int, int, int);
  void FOGSide(double *, int, int, int, int, int);
  void GPR2Side(double *, int, int, int, int, int);
  void Mood(int, int, int, int);

  // Defined in the src/WENO.cpp file
  void WENO(int, int, int, int);

  // Defined in the src/HLL.cpp file
  void HLL();
  void HLLSide(int, int);

  // Defined in the src/Detection.cpp file
  bool Detection(bool);

  void SpaceRecon() {
    for (int quad = 0; quad < nqp; ++quad) {
      for (int var = 0; var < NumVar; ++var) {
        for (int x = 2; x < REdgeX - 2; ++x) {
          for (int y = YStart; y < YEnd; ++y) {
#if SpaceMethod == Weno
            WENO(quad, var, x, y);
#elif SpaceMethod == Fog
            FOG(quad, var, x, y);
#elif SpaceMethod == Gp1
            GPR1(quad, var, x, y);
#elif SpaceMethod == Gp2
            GPR2(quad, var, x, y);
#elif SpaceMethod == Mood53
            Mood(quad, var, x, y);
#endif
          }
        }
>>>>>>> BugResetBranch:2D/HD/include/FluxClass.hpp
      }
    }
  }

<<<<<<< HEAD:include/FluxClass.hpp
  void Quad_Calculation(int Xstart, int Xstop, int Ystart, int Ystop) {
    double R1Vec[5];
    Cell *Center;
    for (int i = Xstart; i < Xstop; ++i) {
      for (int j = Ystart; j < Ystop; ++j) {
        Center = &Cells[i * YDim + j];
        R1Vec[0] = *Center->ENERGY;
        R1Vec[1] = *Center->LCell->ENERGY;
        R1Vec[2] = *Center->TCell->ENERGY;
        R1Vec[3] = *Center->RCell->ENERGY;
        R1Vec[4] = *Center->BCell->ENERGY;

        ENERGY[0][0][j][i] = cblas_ddot(5, R1Vec, 1, SolutionKer.R1_g1R, 1);
        ENERGY[0][1][j][i] = cblas_ddot(5, R1Vec, 1, SolutionKer.R1_g2R, 1);
        ENERGY[0][2][j][i] = cblas_ddot(5, R1Vec, 1, SolutionKer.R1_g1L, 1);
        ENERGY[0][3][j][i] = cblas_ddot(5, R1Vec, 1, SolutionKer.R1_g2L, 1);

        ENERGY[1][0][j][i] = cblas_ddot(5, R1Vec, 1, SolutionKer.R1_g1U, 1);
        ENERGY[1][1][j][i] = cblas_ddot(5, R1Vec, 1, SolutionKer.R1_g2U, 1);
        ENERGY[1][2][j][i] = cblas_ddot(5, R1Vec, 1, SolutionKer.R1_g1D, 1);
        ENERGY[1][3][j][i] = cblas_ddot(5, R1Vec, 1, SolutionKer.R1_g2D, 1);

        R1Vec[0] = *Center->MOMX;
        R1Vec[1] = *Center->LCell->MOMX;
        R1Vec[2] = *Center->TCell->MOMX;
        R1Vec[3] = *Center->RCell->MOMX;
        R1Vec[4] = *Center->BCell->MOMX;

        MOMX[0][0][j][i] = cblas_ddot(5, R1Vec, 1, SolutionKer.R1_g1R, 1);
        MOMX[0][1][j][i] = cblas_ddot(5, R1Vec, 1, SolutionKer.R1_g2R, 1);
        MOMX[0][2][j][i] = cblas_ddot(5, R1Vec, 1, SolutionKer.R1_g1L, 1);
        MOMX[0][3][j][i] = cblas_ddot(5, R1Vec, 1, SolutionKer.R1_g2L, 1);

        MOMX[1][0][j][i] = cblas_ddot(5, R1Vec, 1, SolutionKer.R1_g1U, 1);
        MOMX[1][1][j][i] = cblas_ddot(5, R1Vec, 1, SolutionKer.R1_g2U, 1);
        MOMX[1][2][j][i] = cblas_ddot(5, R1Vec, 1, SolutionKer.R1_g1D, 1);
        MOMX[1][3][j][i] = cblas_ddot(5, R1Vec, 1, SolutionKer.R1_g2D, 1);

        R1Vec[0] = *Center->MOMY;
        R1Vec[1] = *Center->LCell->MOMY;
        R1Vec[2] = *Center->TCell->MOMY;
        R1Vec[3] = *Center->RCell->MOMY;
        R1Vec[4] = *Center->BCell->MOMY;

        MOMY[0][0][j][i] = cblas_ddot(5, R1Vec, 1, SolutionKer.R1_g1R, 1);
        MOMY[0][1][j][i] = cblas_ddot(5, R1Vec, 1, SolutionKer.R1_g2R, 1);
        MOMY[0][2][j][i] = cblas_ddot(5, R1Vec, 1, SolutionKer.R1_g1L, 1);
        MOMY[0][3][j][i] = cblas_ddot(5, R1Vec, 1, SolutionKer.R1_g2L, 1);

        MOMY[1][0][j][i] = cblas_ddot(5, R1Vec, 1, SolutionKer.R1_g1U, 1);
        MOMY[1][1][j][i] = cblas_ddot(5, R1Vec, 1, SolutionKer.R1_g2U, 1);
        MOMY[1][2][j][i] = cblas_ddot(5, R1Vec, 1, SolutionKer.R1_g1D, 1);
        MOMY[1][3][j][i] = cblas_ddot(5, R1Vec, 1, SolutionKer.R1_g2D, 1);

        R1Vec[0] = *Center->DENS;
        R1Vec[1] = *Center->LCell->DENS;
        R1Vec[2] = *Center->TCell->DENS;
        R1Vec[3] = *Center->RCell->DENS;
        R1Vec[4] = *Center->BCell->DENS;

        DENS[0][0][j][i] = cblas_ddot(5, R1Vec, 1, SolutionKer.R1_g1R, 1);
        DENS[0][1][j][i] = cblas_ddot(5, R1Vec, 1, SolutionKer.R1_g2R, 1);
        DENS[0][2][j][i] = cblas_ddot(5, R1Vec, 1, SolutionKer.R1_g1L, 1);
        DENS[0][3][j][i] = cblas_ddot(5, R1Vec, 1, SolutionKer.R1_g2L, 1);

        DENS[1][0][j][i] = cblas_ddot(5, R1Vec, 1, SolutionKer.R1_g1U, 1);
        DENS[1][1][j][i] = cblas_ddot(5, R1Vec, 1, SolutionKer.R1_g2U, 1);
        DENS[1][2][j][i] = cblas_ddot(5, R1Vec, 1, SolutionKer.R1_g1D, 1);
        DENS[1][3][j][i] = cblas_ddot(5, R1Vec, 1, SolutionKer.R1_g2D, 1);
      }
    }
  }
=======
double GetPres(int x, int y) {
  return (GAMMA - 1.0) *
         (Cons[Tidx(Ener, x, y)] -
          0.5 * (Cons[Tidx(MomX, x, y)] * Cons[Tidx(MomX, x, y)] /
                 Cons[Tidx(Dens, x, y)]));
}

double GetPresUin(int x, int y) {
  return (GAMMA - 1.0) * (Uin[Tidx(Ener, x, y)] -
                          0.5 * (Uin[Tidx(MomX, x, y)] * Uin[Tidx(MomX, x, y)] /
                                 Uin[Tidx(Dens, x, y)]));
}

void Recon() {
  for (int quad = 0; quad < nqp; ++quad) {
    for (int var = 0; var < NumVar; ++var) {
      for (int x = XStart; x < XEnd; ++x) {
        for (int y = YStart; y < YEnd; ++y) {
          Cons[Tidx(var, x, y)] -=
              *dt *
              (Flux[quad][var][idx(x, y)] - Flux[quad][var][idx(x - 1, y)]) /
              dx;
        }
      }
    }
  }
}

void ReconSide(int quad, int var, int x, int y) {
  Cons[Tidx(var, x, y)] =
      Uin[Tidx(var, x, y)] -
      *dt * (Flux[quad][var][idx(x, y)] - Flux[quad][var][idx(x - 1, y)]) / dx;
}
}
;

#endif // FLUXCLASS_H_
