#ifndef DOMAINCLASS_H_
#define DOMAINCLASS_H_

#include "GP_Kernel.hpp"
#include "Parameters.h"
#include "definitions.hpp"

#include <algorithm>
#include <cblas.h>
#include <cmath>
#include <iomanip>

class Domain {
public:
  /***********************************************/
  /*********** Internal Data Objects *************/
  /***********************************************/
  double *Dens, *Pres, *Xvel;
  double *DensP, *MomX, *Energy;
  int count;
  double *Cs, *Buffer;
  double *RS_CsL, *RS_CsR;
  double **FluxWalls_Cons;
  double **FluxWalls_Prims;
  double *CellFlux;
  double T, dt, dt_sim;
  double *CopyBuffer;
  double *ConsCopy;
  double *Cons, *Prims;
  double *XDomain;
  bool MoodFinished = true;
  int *MoodOrd;
  int *Troubled;
  int Percentage, CurrentProgress;

  void (Domain::*BC)();
  void (Domain::*IC)();
  void (Domain::*RK_TimeStepper)();
  void (Domain::*SpaceRecon)(int, int);

  GP_Kernel Ker;

  /***********************************************/
  /*************** Class Constructor *************/
  /***********************************************/
  Domain() {

    // Allocate Variables
    XDomain = new double[xDim];
    Cons = new double[NumVar * xDim];
    Prims = new double[NumVar * xDim];
    Buffer = new double[xDim];
    Cs = new double[xDim];
    Percentage = 0;

    FluxWalls_Cons = new double *[2 * NDIMS];
    FluxWalls_Cons[LEFT] = new double[NumVar * xDim];
    FluxWalls_Cons[RIGHT] = new double[NumVar * xDim];

    FluxWalls_Prims = new double *[2];
    FluxWalls_Prims[LEFT] = new double[NumVar * xDim];
    FluxWalls_Prims[RIGHT] = new double[NumVar * xDim];

    CellFlux = new double[NumVar * xDim];
    CopyBuffer = new double[NumVar * xDim];
    Buffer = new double[xDim];
    Cs = new double[xDim];
    RS_CsL = new double[xDim];
    RS_CsR = new double[xDim];

#if SpaceMethod == GPR1
    Ker.calculate_Preds1D(1);
    SpaceRecon = &Domain::GP1;

#elif SpaceMethod == GPR2
    Ker.calculate_Preds1D(2);
    SpaceRecon = &Domain::GP2;

#elif SpaceMethod == MOOD531
    Ker.calculate_Preds1D(1);
    Ker.calculate_Preds1D(2);
    SpaceRecon = &Domain::Mood;

    MoodOrd = new int[xDim];
    Troubled = new bool[xDim];
    ConsCopy = new double[NumVar * xDim];

#elif SpaceMethod == FOG
    SpaceRecon = &Domain::Fog;

#elif SpaceMethod == WENO
    SpaceRecon = &Domain::Weno;
#endif

    Dens = Cons;
    MomX = Dens + xDim;
    Energy = MomX + xDim;

    DensP = Prims;
    Xvel = DensP + xDim;
    Pres = Xvel + xDim;

    T = T0;
    dt_sim = 1E-10;

/**********Member Function Pointers***********/
#if TestProblem == SHUOSHER
    BC = &Domain::ShuOsherBC;
    IC = &Domain::ShuOsherIC;
#elif TestProblem == SHOCKTUBE
    IC = &Domain::ShockTubeIC;
    BC = &Domain::NeumannBC;
#else

#if BCS == NEUMANN
    BCS = &Domain::NeumannBC;
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

  /***********************************************/
  /******* Methods defined in other Files ********/
  /***********************************************/

  // Defined in the VarConvert.cpp file
  void Prims2Cons(double *, double *, int, int);
  int Cons2Prim(double *, double *, int, int);

  void SolvePressure();
  void Press(int x);

  // Defined in the BC.cpp file
  void ShuOsherBC();
  void NeumannBC();

  // Defined in the Find_dt.cpp flie
  void Find_dt();
  void Find_Cs(double *Uin, double *SSVector, int start, int stop);
  // double SRHD_CS(int i);
  // double HD_CS(int i);
  void SignalSpeed(double *Uin, double *SSVector, int i, double &CSL,
                   double &CSR);

  // Defined in the TimeSteppers.cpp file
  void RK3();
  void ForwardEuler();

  // Defined in the IC.cpp file
  void ShuOsherIC();
  void ShockTubeIC();

  // Defined in the IO.cpp file
  void writeResults();

  // Defined in the src/FOG.cpp file
  void Fog(int, int);

  // Defined in the src/GP-FVM.cpp file
  void GP1(int, int);
  void GP2(int, int);
  // void GPR1Side(double *, int, int, int, int);
  // void FOGSide(double *, int, int, int, int);
  // void GPR2Side(double *, int, int, int, int);
  // void Mood(int, int, int);

  // Defined in the src/WENO.cpp file
  void Weno(int, int);

  // Defined in the src/MoodSolve.cpp file
  void Mood(int, int);

  // Defined in the src/HLL.cpp file
  void Hll(int, int);
  void HllSide(int);

  // Defined in the src/UpdateSolution.cpp file
  void Recon(int start, int stop);

  // Defined in the src/Detection.cpp file
  bool Detection(bool);

  // Defined in the src/SR_Flux.cpp file
  void SR_Flux(double *Dest, double *P, double *C, int i, int DestinationIdx);
  void SR_HLL_Flux(double *Dest, double *PrL, double *CL, double *PrR,
                   double *CR, double SL, double SR, int i);

  // Defined in the src/HD_Flux.cpp file
  void HD_Flux(double *Dest, double *P, double *C, int i, int DestinationIdx);
  void HD_HLL_Flux(double *Dest, double *PrL, double *CL, double *PrR,
                   double *CR, double SL, double SR, int i);

  void DomainCopy() { std::copy(Cons, Cons + xDim * NumVar, CopyBuffer); }

  void DomainAdd(double a, double b) {
    cblas_dscal(NumVar * xDim, b, Cons, 1);
    cblas_daxpy(NumVar * xDim, a, CopyBuffer, 1, Cons, 1);
  }

  void PrintProgress() {
    CurrentProgress = (100 * T) / TN;
    if (CurrentProgress > Percentage) {
      std::cout << "We are " << CurrentProgress << "% done, The time is " << T
                << std::endl;
      Percentage = CurrentProgress;
    }
  }
};

#endif // DOMAINCLASS_H_
