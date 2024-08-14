#ifndef DOMAINCLASS_H_
#define DOMAINCLASS_H_

#include "CellClass.hpp"
#include "FluxClass.hpp"
#include "OptionsClass.hpp"
#include <cblas.h>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>
// #include <lapacke.h>

// extern "C" {
// extern void cblas_dscal(int, double, double *, int);
// }

class Domain {
public:
  double *DENS, *PRES, *XVEL, *YVEL, *MOMX, *MOMY, *ENERGY, *Cs, *Buffer;
  double gamma, nx, ny, x0, xN, y0, yN, dx, dy, T, TN, dt, dt_sim, cfl;
  double *CopyBuffer, *CopyCons[4];
  double *XQuads1, *XQuads2, *XQuads3;
  double *YQuads1, *YQuads2, *YQuads3;
  int XStart, YStart, XEnd, YEnd, REdgeX, REdgeY, ngc, ndims;
  int xDim, yDim, ell, nqp;
  Cell *Cells;
  double *Cons[4];
  FluxContainer Flux;
  // double *Prims[4] = {DENS,XVEL,YVEL,PRES};
  // double *Cons[4] = {DENS, XVEL, YVEL, PRES};
  void (Domain::*BC)(std::string);
  void (Domain::*IC)();
  void (Domain::*RK_TimeStepper)();
  bool SlowStart;
  GP_Kernel SolutionKer;

  // Class constructor, Takes in the Opts class to build internal variables
  Domain(opt Opts) {
    xDim = Opts.ngc * 2 + Opts.nx;
    yDim = Opts.ngc * 2 + Opts.ny;
    DENS = new double[xDim * yDim];
    PRES = new double[xDim * yDim];
    XVEL = new double[xDim * yDim];
    YVEL = new double[xDim * yDim];
    MOMX = new double[xDim * yDim];
    MOMY = new double[xDim * yDim];
    ENERGY = new double[xDim * yDim];
    CopyBuffer = new double[xDim * yDim * 4];

    CopyCons[0] = CopyBuffer;
    CopyCons[1] = &CopyBuffer[xDim * yDim];
    CopyCons[2] = &CopyBuffer[xDim * yDim * 2];
    CopyCons[3] = &CopyBuffer[xDim * yDim * 3];

    // double *Prims[4] = {DENS, XVEL, YVEL, PRES};
    Cons[0] = DENS;
    Cons[1] = MOMX;
    Cons[2] = MOMY;
    Cons[3] = ENERGY;
    Buffer = new double[xDim * yDim];
    Cs = new double[xDim * yDim];
    Cells = new Cell[xDim * yDim];
    XStart = Opts.ngc;
    YStart = Opts.ngc;
    REdgeX = xDim;
    REdgeY = yDim;
    XEnd = REdgeX - Opts.ngc;
    YEnd = REdgeY - Opts.ngc;
    ngc = Opts.ngc;
    x0 = Opts.x0;
    y0 = Opts.y0;
    xN = Opts.xN;
    yN = Opts.yN;
    nx = Opts.nx;
    ny = Opts.ny;
    ell = Opts.ell;
    dx = (xN - x0) / nx;
    dy = (yN - y0) / ny;
    T = Opts.T0;
    dt_sim = 1E-10;
    TN = Opts.TN;
    nqp = Opts.ngp;
    cfl = Opts.CFL;
    ndims = Opts.Ndims;
    SlowStart = Opts.SlowStart;
    Flux.FluxContainer_init(nqp, REdgeX, REdgeY, ell, dx, dy, Cells);
    SolutionKer.GP_Kernel_init(ell, dx, dy);
  }

  void AssignCells() {
    for (int i = 0; i < REdgeX; ++i) {
      for (int j = 0; j < REdgeY; ++j) {
        // Each cell stores memory locations corresponding to its position in
        // this array.
        // This approach, while memory intensive, should prevent
        // copying memory.

        int Tmp = i * REdgeY + j;
        Cells[Tmp].SetUpCell(i, j, DENS + Tmp, PRES + Tmp, XVEL + Tmp,
                             YVEL + Tmp, MOMX + Tmp, MOMY + Tmp, ENERGY + Tmp,
                             Cs + Tmp, dx, dy, &dt, gamma, &SolutionKer, nqp);

        // std::cout << Cells[Tmp].y << std::endl;

        // Cells[i * REdgeY + j].DENS = &DENS[i * REdgeY + j];
        // Cells[i * REdgeY + j].PRES = &PRES[i * REdgeY + j];
        // Cells[i * REdgeY + j].XVEL = &XVEL[i * REdgeY + j];
        // Cells[i * REdgeY + j].YVEL = &YVEL[i * REdgeY + j];
        // Cells[i * REdgeY + j].MOMX = &MOMX[i * REdgeY + j];
        // Cells[i * REdgeY + j].MOMY = &MOMY[i * REdgeY + j];
        // Cells[i * REdgeY + j].ENERGY = &ENERGY[i * REdgeY + j];
        // Cells[i * REdgeY + j].Cs = &Cs[i * REdgeY + j];
        // Cells[i * REdgeY + j].dx = dx;
        // Cells[i * REdgeY + j].dy = dy;
        // Cells[i * REdgeY + j].dt = &dt;
        // Cells[i * REdgeY + j].gamma = gamma;
        // Cells[i * REdgeY + j].GP = &SolutionKer;
        // Cells[i * REdgeY + j].nqp = nqp;
        // Cells[i * REdgeY + j].x = i;
        // Cells[i * REdgeY + j].y = j;

        // Access points to get left and right cells as well.
        if (i > 0) {
          Cells[i * REdgeY + j].LCell = GetCell(i - 1, j);
        }
        if (j > 0) {
          Cells[i * REdgeY + j].BCell = GetCell(i, j - 1);
        }

        if (i < xDim - 1) {
          Cells[i * REdgeY + j].RCell = GetCell(i + 1, j);
        }
        if (j < REdgeY - 1) {
          Cells[i * REdgeY + j].TCell = GetCell(i, j + 1);
        }
      }
    }

    for (int i = 2; i < REdgeX - 2; ++i) {
      for (int j = 2; j < REdgeY - 2; ++j) {
        (*GetCell(i, j)).Assign_Stencils();
      }
    }
  }

  Cell *GetCell(int x, int y) { return &Cells[x * REdgeY + y]; }

  // Defined in the VarConvert.cpp file
  void Prims2Cons();
  void Cons2Prim();
  void SolvePressure();

  // Defined in the BC.cpp file
  void NeumannBC(std::string);
  void ShuOsherBC(std::string);

  // Defined in the IC.cpp file
  void ShuOsherIC();

  // Defined in the Find_dt.cpp flie
  void Find_dt();
  void Find_Cs();

  // Defined in the TimeSteppers.cpp file
  void RK3();
  void ForwardEuler();

  // Defined in Flux.cpp file
  void Calculate_Quad_Points();
  void SolveRiemann();
  void TimeStep();

  // Definied in the IO.cpp file
  void WriteOutSolution();

  void SaveDomain() {
    for (int i = 0; i < 4; ++i) {
      cblas_dcopy(xDim * yDim, Cons[i], 1, CopyCons[i], 1);
      // CopyBuffer[i] = DENS[i];
      // std::cout << &DENS[i] << std::endl;
    }

    // std::cout << Cons[1] << std::endl;
    // std::cout << DENS << std::endl;
    // std::cout << &Cons[1] << MOMX << std::endl;
  }

  void DomainAdd(double a, double b) {
    // Vector add A*U + B*U^n+1/2;
    for (int i = 0; i < 4; ++i) {
      cblas_dscal(xDim * yDim, a, Cons[i], 1);
    }

    for (int i = 0; i < 4; ++i) {
      cblas_daxpy(xDim * yDim, b, CopyCons[i], 1, Cons[i], 1);
    }
  }

  void Check() {
    for (int i = 0; i < REdgeX; ++i) {
      std::cout << std::setprecision(16) << std::fixed << DENS[i * 7 + 3]
                << std::endl;
    }
    exit(0);
  }

  void TestLapacke(double val) { cblas_dscal(25, val, DENS, 1); }
};

#endif // DOMAINCLASS_H_
