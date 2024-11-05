#ifndef DOMAINCLASS_H_
#define DOMAINCLASS_H_

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
  double *Dens, *DensP, *Pres, *Xvel, *Yvel, *Zvel;
  double *MomX, *MomY, *MomZ, *Energy, *Cs, *Buffer;
  double *RS_CsL, *RS_CsR;
  double **FluxWalls_Cons;
  double **FluxWalls_Prims;
  double *CellFlux;
  double T, dt, dt_sim;
  double *CopyBuffer;
  double *ConsCopy;
  double *U1, *U2, *U3, *U4, *Fl, *FU3, *FU4;
  double *Cons, *Prims;
  bool MoodFinished = true;
  int *MoodOrd;
  int *Troubled;

  void (Domain::*BC)();
  void (Domain::*IC)();
  void (Domain::*RiemannSolver)(int, int);
  void (Domain::*RK_TimeStepper)();
  void (Domain::*SpaceRecon)(int, int);

  GP_Kernel Ker;

  /***********************************************/
  /*************** Class Constructor *************/
  /***********************************************/
  Domain() {

    // Allocate Variables
    Cons = new double[NumVar * xDim];
    Prims = new double[NumVar * xDim];
    Buffer = new double[xDim];
    Cs = new double[xDim];

    FluxWalls_Cons = new double *[2];
    FluxWalls_Cons[LEFT] = new double[NumVar * xDim];
    FluxWalls_Cons[RIGHT] = new double[NumVar * xDim];

    FluxWalls_Prims = new double *[2];
    FluxWalls_Prims[LEFT] = new double[NumVar * xDim];
    FluxWalls_Prims[RIGHT] = new double[NumVar * xDim];

    CellFlux = new double[NumVar * xDim];
    CopyBuffer = new double[NumVar * xDim];

#if RK_Method > 3
    U1 = new double[NumVar * xDim];
    U2 = new double[NumVar * xDim];
    U3 = new double[NumVar * xDim];
    U4 = new double[NumVar * xDim];
    FU3 = new double[NumVar * xDim];
    FU4 = new double[NumVar * xDim];
    Fl = new double[NumVar * xDim];
#endif

    Buffer = new double[xDim];
    Cs = new double[xDim];
    RS_CsL = new double[xDim];
    RS_CsR = new double[xDim];

#if SpaceMethod == GPR1
    SolutionKer.calculate_Preds1D(1);
    SpaceRecon = &Domain::GP1;

#elif SpaceMethod == GPR2
    SolutionKer.calculate_Preds1D(2);
    SpaceRecon = &Domain::GP2;

#elif SpaceMethod == MOOD531
    SolutionKer.calculate_Preds1D(1);
    SolutionKer.calculate_Preds1D(2);
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
    MomY = MomX + xDim;
    MomZ = MomY + xDim;
    Energy = MomZ + xDim;

    DensP = Prims;
    Xvel = DensP + xDim;
    Yvel = Xvel + xDim;
    Zvel = Yvel + xDim;
    Pres = Zvel + xDim;

    T = T0;
    dt_sim = 1E-10;

/**********Member Function Pointers***********/
#if TestProblem == SHUOSHER
    BC = &Domain::ShuOsherBC;
    IC = &Domain::ShuOsherIC;
#elif TestProblem == SHOCKTUBE
    IC = &Domain::ShockTubeIC;
    BC = &Domain::ShuOsherBC;
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
#elif RK_Method == 4
    RK_TimeStepper = &Domain::RK4;
#endif

#if RIEMANN == HLL
    RiemannSolver = &Domain::Hll;
#elif RIEMANN == HLLC
    RiemannSolver = &Domain::Hllc;
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
  void RK4();

  // Defined in the IC.cpp file
  void ShuOsherIC();
  void ShockTubeIC();

  // Defined in the EnergyInverter.cpp file
  int EnergyInverter(double *Uin, double *Uout, int i);

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

  // Defined in the src/HLLC.cpp file
  void Hllc(int Start, int Stop);

  // Defined in the src/UpdateSolution.cpp file
  void Recon(int start, int stop);

  // Defined in the src/Detection.cpp file
  bool Detection(bool);

  // Defined in the src/NewtonPressureFinder.cpp file
  int NaiveNewton(double *Uin, double *Uout, int i);

  // Defined in the src/PressureFix.cpp file
  int PressureFix(double *Uin, double *Uout, int i);

  // Defined in the src/SR_Flux.cpp file
  void SR_Flux(double *Dest, double *P, double *C, int i, int DestinationIdx);
  void SR_HLL_Flux(double *Dest, double *PrL, double *CL, double *PrR,
                   double *CR, double SL, double SR, int i);

  void DomainCopy(double *Uin, double *Ufill) {
    std::copy(Uin, Uin + xDim * NumVar, Ufill);
  }

  void DomainAdd(double a, double b, double *Uin, double *Uresult) {
    cblas_dscal(NumVar * xDim, b, Uresult, 1);
    cblas_daxpy(NumVar * xDim, a, Uin, 1, Uresult, 1);
  }
};

#endif // DOMAINCLASS_H_
