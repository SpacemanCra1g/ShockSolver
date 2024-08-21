#ifndef DOMAINCLASS_H_
#define DOMAINCLASS_H_

#include "CellClass.hpp"
#include "FluxClass.hpp"
#include "GP_Kernel.hpp"
#include "Parameters.h"

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
  /***********************************************/
  /*********** Internal Data Objects *************/
  /***********************************************/
  double *DENS, *PRES, *XVEL, *YVEL, *MOMX, *MOMY, *ENERGY, *Cs, *Buffer;
  double T, dt, dt_sim;
  double *CopyBuffer, *CopyCons;

  Cell *Cells;
  double *Cons;
  double *Prims[NumVar];
  FluxClass Flux;
  // double *Prims[4] = {DENS,XVEL,YVEL,PRES};
  void (Domain::*BC)(std::string);
  void (Domain::*IC)();
  void (Domain::*RK_TimeStepper)();

  GP_Kernel SolutionKer;

  /***********************************************/
  /******* Methods defined in other Files ********/
  /***********************************************/

  // Defined in the VarConvert.cpp file
  void Prims2Cons();
  void Cons2Prim();
  void SolvePressure();

  // Defined in the BC.cpp file
  void ShuOsherBC(std::string);
  void NeumannBC(std::string);

  // Defined in the Find_dt.cpp flie
  void Find_dt();
  void Find_Cs();

  // Defined in the TimeSteppers.cpp file
  void RK3();
  void ForwardEuler();

  // Defined in the IC.cpp file
  void ShuOsherIC();

  /***********************************************/
  /*************** Class Constructor *************/
  /***********************************************/
  Domain() {

    Buffer = new double[xDim * yDim];
    Cs = new double[xDim * yDim];
    Cells = new Cell[xDim * yDim];

    T = T0;
    dt_sim = 1E-10;

    SolutionKer.GP_Kernel_init();
    Flux.Fluxinit(Cons);

    Cons = new double[xDim * yDim * NumVar];
    DENS = Cons;
    PRES = new double[xDim * yDim];
    XVEL = new double[xDim * yDim];
    YVEL = new double[xDim * yDim];
    MOMX = Cons + xDim * yDim;
    ENERGY = MOMX + yDim * xDim;

#if NumVar > 3
    MOMY = ENERGY + yDim * xDim;
#endif

    CopyBuffer = new double[xDim * yDim * 4];

    double *CopyCons[4] = {CopyBuffer, &CopyBuffer[xDim * yDim],
                           &CopyBuffer[xDim * yDim * 2],
                           &CopyBuffer[xDim * yDim * 3]};

    double *Prims[4] = {DENS, XVEL, YVEL, PRES};

/**********Member Function Pointers***********/
#if TestProblem == ShuOsher
    BC = &Domain::ShuOsherBC;
    IC = &Domain::ShuOsherIC;
#else

#if BC == Neumann
    BC = &Domain::NeumannBC;
#else
    std::cout << "Invalid Boundary Conditions \nExiting" << std::endl;
    exit(0);
#endif
#endif

#if RK_Method == 1
    RK_TimeStepper = &Domain::ForwardEuler;
#elif RK_Method == 3
    RK_TimeStepper = &Domain::RK3;
#endif
  }
  void AssignCells() {
    for (int i = 0; i < xDim; ++i) {
      for (int j = 0; j < yDim; ++j) {
        // Each cell stores memory locations corresponding to its position in
        // this array.
        // This approach, while memory intensive, should prevent
        // copying memory.
        Cells[idx(i, j)].DENS = &DENS[idx(i, j)];
        Cells[idx(i, j)].PRES = &PRES[idx(i, j)];
        Cells[idx(i, j)].XVEL = &XVEL[idx(i, j)];
        Cells[idx(i, j)].YVEL = &YVEL[idx(i, j)];
        Cells[idx(i, j)].MOMX = &MOMX[idx(i, j)];
        Cells[idx(i, j)].MOMY = &MOMY[idx(i, j)];
        Cells[idx(i, j)].ENERGY = &ENERGY[idx(i, j)];
        Cells[idx(i, j)].Cs = &Cs[idx(i, j)];
        Cells[idx(i, j)].x = i;
        Cells[idx(i, j)].y = j;
        Cells[idx(i, j)].GP_Weight = &SolutionKer;

        // Access points to get left and right cells as well.
        if (i > 0) {
          Cells[idx(i, j)].LCell = &Cells[idx(i - 1, j)];
        }
        if (j > 0) {
          Cells[idx(i, j)].BCell = &Cells[idx(i, j - 1)];
        }

        if (i > xDim - 1) {
          Cells[idx(i, j)].RCell = &Cells[idx(i + 1, j)];
        }
        if (j > yDim - 1) {
          Cells[idx(i, j)].TCell = &Cells[idx(i, j + 1)];
        }
      }
    }
  }

  void DomainCopy() {
    std::copy(Cons, Cons + yDim * xDim * NumVar, CopyBuffer);
  }

  void DomainAdd(double a, double b) {
    cblas_dscal(NumVar * xDim * yDim, b, Cons, 1);
    cblas_daxpy(NumVar * xDim * yDim, a, CopyBuffer, 1, Cons, 1);
  }

  Cell *GetCell(int x, int y) { return &Cells[idx(x, y)]; }

  void TestLapacke(double val) { cblas_dscal(25, val, DENS, 1); }
};

#endif // DOMAINCLASS_H_
