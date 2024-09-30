#ifndef DOMAINCLASS_H_
#define DOMAINCLASS_H_

#include "FluxClass.hpp"
#include "GP_Kernel.hpp"
#include "Parameters.h"

#include <algorithm>
#include <cblas.h>
#include <cmath>

class Domain {
public:
  /***********************************************/
  /*********** Internal Data Objects *************/
  /***********************************************/
  double *DENS, *DENSP, *PRES, *XVEL, *YVEL, *ZVEL;
  double *MOMX, *MOMY, *MOMZ, *ENERGY, *Cs, *Buffer;
  double T, dt, dt_sim;
  double *CopyBuffer;
  double *ConsCopy;
  bool MoodFinished = true;

  double *Cons;
  double *Prims;

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
  void Press(int x);

  // Defined in the BC.cpp file
  void ShuOsherBC(std::string);
  void NeumannBC(std::string);

  // Defined in the Find_dt.cpp flie
  void Find_dt();
  void Find_Cs();
  double SRHD_CS(int i);
  double HD_CS(int i);

  // Defined in the TimeSteppers.cpp file
  void RK3();
  void ForwardEuler();

  // Defined in the IC.cpp file
  void ShuOsherIC();
  void ShockTubeIC();

  // Defined in the IO.cpp file
  void writeResults();

  /***********************************************/
  /*************** Class Constructor *************/
  /***********************************************/
  Domain() {

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
    Prims = new double[xDim * NumVar];

    DENS = Cons;
    MOMX = Cons + xDim;
    MOMY = Cons + 2 * xDim;
    MOMZ = Cons + 3 * xDim;
    ENERGY = Cons + 4 * xDim;

    DENSP = Prims;
    XVEL = DENSP + xDim;
    YVEL = XVEL + xDim;
    ZVEL = YVEL + xDim;
    PRES = ZVEL + xDim;

    CopyBuffer = new double[xDim * NumVar];
    Flux.Fluxinit(Cons, &dt);

#if SpaceMethod == Mood53
    ConsCopy = new double[xDim * NumVar];
    Flux.Uin = ConsCopy;
#endif

/**********Member Function Pointers***********/
#if TestProblem == SHUOSHER
    BC = &Domain::ShuOsherBC;
    IC = &Domain::ShuOsherIC;
#elif TestProblem == SHOCKTUBE
    IC = &Domain::ShockTubeIC;
    BC = &Domain::ShuOsherBC;
#else

#if BCs == NEUMANN
    BCs = &Domain::NeumannBC;
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

  void DomainCopy() { std::copy(Cons, Cons + xDim * NumVar, CopyBuffer); }

  void DomainAdd(double a, double b) {
    cblas_dscal(NumVar * xDim, b, Cons, 1);
    cblas_daxpy(NumVar * xDim, a, CopyBuffer, 1, Cons, 1);
  }
};

#endif // DOMAINCLASS_H_
