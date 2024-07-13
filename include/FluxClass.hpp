#ifndef FLUXCLASS_H_
#define FLUXCLASS_H_

#include "CellClass.hpp"
#include "GP_Kernel.hpp"
#include <vector>

class FluxContainer {

public:
  std::vector<double> **ENERGY[2];
  std::vector<double> **MOMX[2];
  std::vector<double> **MOMY[2];
  std::vector<double> **DENS[2];
  int XDim;
  int YDim;
  GP_Kernel SolutionKer;
  Cell *Cells;

  void FluxContainer_init(const int nqp, int REdgeX, int REdgeY, int ell,
                          double dx, double dy, Cell *C) {

    ENERGY[0] = new std::vector<double> *[2 * nqp];
    ENERGY[1] = new std::vector<double> *[2 * nqp];

    MOMX[0] = new std::vector<double> *[2 * nqp];
    MOMX[1] = new std::vector<double> *[2 * nqp];

    MOMY[0] = new std::vector<double> *[2 * nqp];
    MOMY[1] = new std::vector<double> *[2 * nqp];

    DENS[0] = new std::vector<double> *[2 * nqp];
    DENS[1] = new std::vector<double> *[2 * nqp];

    for (int i = 0; i < 2 * nqp; ++i) {
      ENERGY[0][i] = new std::vector<double>[REdgeY];
      ENERGY[1][i] = new std::vector<double>[REdgeX];

      MOMX[0][i] = new std::vector<double>[REdgeY];
      MOMX[1][i] = new std::vector<double>[REdgeX];

      MOMY[0][i] = new std::vector<double>[REdgeY];
      MOMY[1][i] = new std::vector<double>[REdgeX];

      DENS[0][i] = new std::vector<double>[REdgeY];
      DENS[1][i] = new std::vector<double>[REdgeX];
    }

    Cells = C;

    XDim = REdgeX;
    YDim = REdgeY;

    SolutionKer.GP_Kernel_init(ell, dx, dy);

    for (int j = 0; j < 2 * nqp; ++j) {

      for (int i = 0; i < REdgeX; ++i) {
        ENERGY[1][j][i].resize(REdgeY - 1);
        MOMX[1][j][i].resize(REdgeY - 1);
        MOMY[1][j][i].resize(REdgeY - 1);
        DENS[1][j][i].resize(REdgeY - 1);
      }

      for (int i = 0; i < REdgeY; ++i) {
        ENERGY[0][j][i].resize(REdgeX - 1);
        MOMX[0][j][i].resize(REdgeX - 1);
        MOMY[0][j][i].resize(REdgeX - 1);
        DENS[0][j][i].resize(REdgeX - 1);
      }
    }
  }

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
};

#endif // FLUXCLASS_H_
