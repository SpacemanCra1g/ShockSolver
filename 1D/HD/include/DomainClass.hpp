#ifndef DOMAINCLASS_H_
#define DOMAINCLASS_H_

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
  double *CopyBuffer;
  double *ConsCopy;
  int count;
  bool MoodFinished = true;

  double *Cons;
  double *Prims[NumVar];

  FluxClass Flux;

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
  void SodIC();

  // Defined in the IO.cpp file
  void writeResults();

  /***********************************************/
  /*************** Class Constructor *************/
  /***********************************************/
  Domain() {

    // Prims[Dens] = DENS;
    // Prims[VelX] = XVEL;
    // Prims[Pres] = PRES;

#if SpaceMethod == Gp1
    SolutionKer.calculate_Preds1D(1);
    Flux.Kern = &SolutionKer;
#elif SpaceMethod == Gp2
    SolutionKer.calculate_Preds1D(2);
    Flux.Kern = &SolutionKer;
#elif SpaceMethod == Mood53
    SolutionKer.calculate_Preds1D(1);
    SolutionKer.calculate_Preds1D(2);
    Flux.Kern = &SolutionKer;
    Flux.MoodOrd = new int[xDim];
    Flux.Troubled = new bool[xDim];
    std::fill(Flux.MoodOrd, Flux.MoodOrd + xDim, 5);
    std::fill(Flux.Troubled, Flux.Troubled + xDim, false);
#endif

    Buffer = new double[xDim];
    Cs = new double[xDim];

    T = T0;
    dt_sim = 1E-10;

    // SolutionKer.GP_Kernel_init();

    Cons = new double[xDim * NumVar];
    DENS = &Cons[0];
    PRES = new double[xDim];
    XVEL = new double[xDim];

    MOMX = &Cons[MomX * xDim];
    ENERGY = &Cons[Ener * (xDim)];

    CopyBuffer = new double[xDim * NumVar];
    Flux.Fluxinit(Cons, &dt);

#if SpaceMethod == Mood53
    ConsCopy = new double[xDim * NumVar];
    Flux.Uin = ConsCopy;
#endif

/**********Member Function Pointers***********/
#if TestProblem == ShuOsher
    BC = &Domain::ShuOsherBC;
    IC = &Domain::ShuOsherIC;
#elif TestProblem == Sod
    BC = &Domain::NeumannBC;
    IC = &Domain::SodIC;
#else

#endif

#if RK_Method == 1
    RK_TimeStepper = &Domain::ForwardEuler;
#elif RK_Method == 3
    RK_TimeStepper = &Domain::RK3;
#endif
  }

  void DomainCopy() { std::copy(Cons, Cons + xDim * NumVar, CopyBuffer); }

  void DomainAdd(double a, double b) {
    cblas_dscal(NumVar * xDim, b, Cons, 1);
    cblas_daxpy(NumVar * xDim, a, CopyBuffer, 1, Cons, 1);
  }
};

#endif // DOMAINCLASS_H_
