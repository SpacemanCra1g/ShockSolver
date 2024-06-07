#ifndef DOMAINCLASS_H_
#define DOMAINCLASS_H_

#include "CellClass.hpp"
#include "OptionsClass.hpp"
#include <cblas.h>
#include <cmath>
#include <iostream>
// #include <lapacke.h>

// extern "C" {
// extern void cblas_dscal(int, double, double *, int);
// }

class Domain {
public:
  double *DENS, *PRES, *XVEL, *YVEL, *MOMX, *MOMY, *ENERGY, *Cs, *Buffer;
  double gamma, nx, ny, x0, xN, y0, yN, dx, dy, T, TN, dt, dt_sim, cfl;
    double *CopyBuffer;
  int XStart, YStart, XEnd, YEnd, REdgeX, REdgeY, ngc, ndims;
  int xDim, yDim;
  Cell *Cells;
    // double *Prims[4] = {DENS,XVEL,YVEL,PRES};
  void (Domain::*BC)(std::string);
  void (Domain::*IC)();
  void (Domain::*RK_TimeStepper)();
  bool SlowStart;

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
    CopyBuffer = new double[xDim*yDim*4];
    double *CopyCons[4] ={CopyBuffer,&CopyBuffer[xDim*yDim],
                   &CopyBuffer[xDim*yDim*2],&CopyBuffer[xDim*yDim*3]};
    double *Prims[4] = {DENS,XVEL,YVEL,PRES};
    double *Cons[4] = {DENS,MOMX,MOMY,ENERGY};
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
    dx = (xN - x0) / nx;
    dy = (yN - y0) / ny;
    T = Opts.T0;
    dt_sim = 1E-10;
    TN = Opts.TN;
    cfl = Opts.CFL;
    ndims = Opts.Ndims;
    SlowStart = Opts.SlowStart;
  }
  void AssignCells() {
    for (int i = 0; i < xDim; ++i) {
      for (int j = 0; j < yDim; ++j) {
        // Each cell stores memory locations corresponding to it's position in
        // this array This approach, while memory intensive, should prevent
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
