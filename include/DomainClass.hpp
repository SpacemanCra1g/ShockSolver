#ifndef DOMAINCLASS_H_
#define DOMAINCLASS_H_

#include "OptionsClass.hpp"
#include <cblas.h>
#include <iostream>
// #include <lapacke.h>

// extern "C" {
// extern void cblas_dscal(int, double, double *, int);
// }

class Cell {
public:
  double *DENS, *PRES, *XVEL, *YVEL, *MOMX, *MOMY, *ENERGY, *gamma;
  Cell *LCell, *RCell, *TCell, *BCell;
  int x, y;

  // Defined VarConvert.cpp
  void Cons2Prims();
  void Prims2Cons();
  double GetPres();

  void Assign(std::string Var, double val) {
    if (Var == "DENS") {
      *DENS = val;
    } else if (Var == "PRES") {
      *PRES = val;
    } else if (Var == "XVEL") {
      *XVEL = val;
    } else if (Var == "YVEL") {
      *YVEL = val;
    }
  }
};

class Domain {
public:
  double *DENS, *PRES, *XVEL, *YVEL, *MOMX, *MOMY, *ENERGY;
  double gamma, nx, ny, x0, xN, y0, yN, dx, dy;
  int XStart, YStart, XEnd, YEnd, REdgeX, REdgeY, ngc;

public:
  int xDim, yDim;
  Cell *Cells;
  void (Domain::*BC)(std::string);
  void (Domain::*IC)();

  Domain(opt &Opts) {
    xDim = Opts.ngc * 2 + Opts.nx;
    yDim = Opts.ngc * 2 + Opts.ny;
    DENS = new double[xDim * yDim];
    PRES = new double[xDim * yDim];
    XVEL = new double[xDim * yDim];
    YVEL = new double[xDim * yDim];
    MOMX = new double[xDim * yDim];
    MOMY = new double[xDim * yDim];
    ENERGY = new double[xDim * yDim];
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
    dx = (x0 - xN) / nx;
    dy = (y0 - xN) / ny;
  }
  void AssignCells() {
    for (int i = 0; i < xDim; ++i) {
      for (int j = 0; j < yDim; ++j) {
        Cells[i * yDim + j].DENS = &DENS[i * yDim + j];
        Cells[i * yDim + j].PRES = &PRES[i * yDim + j];
        Cells[i * yDim + j].XVEL = &XVEL[i * yDim + j];
        Cells[i * yDim + j].YVEL = &YVEL[i * yDim + j];
        Cells[i * yDim + j].MOMX = &MOMX[i * yDim + j];
        Cells[i * yDim + j].MOMY = &MOMY[i * yDim + j];
        Cells[i * yDim + j].ENERGY = &ENERGY[i * yDim + j];
        Cells[i * yDim + j].gamma = &gamma;
        Cells[i * yDim + j].x = i;
        Cells[i * yDim + j].y = j;

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

  // Defined in the BC.cpp file
  void NeumannBC(std::string);
  void ShuOsherBC(std::string);

  // Defined in the IC.cpp file
  void ShuOsherIC();

  void TestLapacke(double val) { cblas_dscal(25, val, DENS, 1); }
};

#endif // DOMAINCLASS_H_
