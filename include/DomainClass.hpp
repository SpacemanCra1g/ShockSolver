#ifndef DOMAINCLASS_H_
#define DOMAINCLASS_H_

#include "CellClass.hpp"
#include "FluxClass.hpp"
#include "GP_Kernel.hpp"
#include "OptionsClass.hpp"
#include <algorithm>
#include <cblas.h>
#include <cmath>
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
  double *CopyBuffer, *CopyCons;
  double *XQuads1, *XQuads2, *XQuads3;
  double *YQuads1, *YQuads2, *YQuads3;
  int XStart, YStart, XEnd, YEnd, REdgeX, REdgeY, ngc, ndims;
  int xDim, yDim, ell, nqp;
  Cell *Cells;
  double *Cons;
  FluxClass Flux;
  // double *Prims[4] = {DENS,XVEL,YVEL,PRES};
  void (Domain::*BC)(std::string);
  void (Domain::*IC)();
  void (Domain::*RK_TimeStepper)();
  bool SlowStart;
  bool twoD;
  GP_Kernel SolutionKer;

  // Class constructor, Takes in the Opts class to build internal variables
  Domain(opt Opts) {
    if (Opts.Ndims > 1) {
      twoD = true;
    }
    xDim = Opts.ngc * 2 + Opts.nx;
    if (twoD) {
      yDim = Opts.ngc * 2 + Opts.ny;
    } else {
      yDim = 1;
    }

    Cons = new double[xDim * yDim * 4];
    DENS = Cons;
    PRES = new double[xDim * yDim];
    XVEL = new double[xDim * yDim];
    YVEL = new double[xDim * yDim];
    MOMX = Cons + xDim * yDim;
    MOMY = MOMX + yDim * xDim;
    ENERGY = MOMY + yDim * xDim;

    CopyBuffer = new double[xDim * yDim * 4];

    double *CopyCons[4] = {CopyBuffer, &CopyBuffer[xDim * yDim],
                           &CopyBuffer[xDim * yDim * 2],
                           &CopyBuffer[xDim * yDim * 3]};

    double *Prims[4] = {DENS, XVEL, YVEL, PRES};

    // Not sure what these are for, remember to look into them
    XQuads1 = new double[(xDim - 1) * yDim];
    YQuads1 = new double[xDim * (yDim - 1)];

    Buffer = new double[xDim * yDim];
    Cs = new double[xDim * yDim];
    Cells = new Cell[xDim * yDim];
    XStart = Opts.ngc;
    if (twoD) {
      YStart = Opts.ngc;
      REdgeY = yDim;
      YEnd = REdgeY - Opts.ngc;
    } else {
      YStart = 0;
      REdgeY = 1;
      YEnd = 1;
    }
    REdgeX = xDim;

    XEnd = REdgeX - Opts.ngc;

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
    SolutionKer.GP_Kernel_init(ell, dx, dy);
    Flux.Fluxinit(ndims, XStart, XEnd, YStart, YEnd, nqp, Cons);
  }
  void AssignCells() {
    for (int i = 0; i < xDim; ++i) {
      for (int j = 0; j < yDim; ++j) {
        // Each cell stores memory locations corresponding to its position in
        // this array.
        // This approach, while memory intensive, should prevent
        // copying memory.
        Cells[i * yDim + j].DENS = &DENS[i * yDim + j];
        Cells[i * yDim + j].PRES = &PRES[i * yDim + j];
        Cells[i * yDim + j].XVEL = &XVEL[i * yDim + j];
        Cells[i * yDim + j].YVEL = &YVEL[i * yDim + j];
        Cells[i * yDim + j].MOMX = &MOMX[i * yDim + j];
        Cells[i * yDim + j].MOMY = &MOMY[i * yDim + j];
        Cells[i * yDim + j].ENERGY = &ENERGY[i * yDim + j];
        Cells[i * yDim + j].Cs = &Cs[i * yDim + j];
        Cells[i * yDim + j].gamma = &gamma;
        Cells[i * yDim + j].x = i;
        Cells[i * yDim + j].y = j;
        Cells[i * yDim + j].GP_Weight = &SolutionKer;

        // Access points to get left and right cells as well.
        if (i > 0) {
          Cells[i * yDim + j].LCell = &Cells[(i - 1) * yDim + j];
        }
        if (j > 0) {
          Cells[i * yDim + j].BCell = &Cells[(i)*yDim + j - 1];
        }

        if (i > xDim - 1) {
          Cells[i * yDim + j].RCell = &Cells[(i + 1) * yDim + j];
        }
        if (j > yDim - 1) {
          Cells[i * yDim + j].TCell = &Cells[(i)*yDim + j + 1];
        }
      }
    }
  }

  void DomainCopy() { std::copy(Cons, Cons + yDim * xDim * 4, CopyCons); }
  void DomainAdd(double a, double b) {
    cblas_dscal(4 * xDim * yDim, b, Cons, 1);
    cblas_daxpy(4 * xDim * yDim, a, CopyCons, 1, Cons, 1);
  }

  Cell *GetCell(int x, int y) { return &Cells[x * yDim + y]; }

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

  void TestLapacke(double val) { cblas_dscal(25, val, DENS, 1); }
};

#endif // DOMAINCLASS_H_
