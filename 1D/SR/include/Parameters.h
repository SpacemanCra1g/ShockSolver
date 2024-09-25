#include "definitions.hpp"

#define NDIMS 1
#define NGC 3

/****************
 * X Parameters *
 ****************/
#define NX 50
#define X0 0.0
#define XN 1.0

#define RhoL 10.0
#define RhoR 1.0

#define PL (40.0 / 3.0)
#define PR (0.000001 * (2.0 / 3.0))

#define XVelL 0.0
#define YVelL 0.0
#define ZVelL 0.0

#define XVelR 0.0
#define YVelR 0.0
#define ZVelR 0.0

/****************
 * Y Parameters *
 ****************/

#define NumVar 5

/***************
 * T Parameters *
 ****************/
#define T0 0.0
#define TN 1.8

/*****************
 *Run Parameters *
 *****************/
#define SpaceMethod Weno

#define GAMMA (5.0 / 3.0)
#define CFL 0.8
#define TestProblem SHOCKTUBE
#define BCs NEUMANN
#define EOS IdealGas

#define RK_Method 1
#define ell 6.0
#define MoodOrder 5
#define SlowStart True

#define nqp 1

/****************
 *NN Parameters *
 ****************/
#define Use_NN False
#define UseDMP True
#define NegativeSlope 0.01
#define NN_Thresh .99

/******************************
 * ****************************
 * ****************************
 * Compiler defined constants *
 * Don't touch                *
 ******************************/
#define xDim (2 * NGC + NX)
#define XStart NGC
#define REdgeX xDim
#define XEnd (REdgeX - NGC)
#define dx ((XN - X0) / NX)

#define sigdel (ell * 1.4142135623730950488016887240L)

#define Dens 0
#define MomX 1
#define MomY 2
#define MomZ 3
#define Ener 4

#define DensP 0
#define VelX 1
#define VelY 2
#define VelZ 3
#define Pres 4

#define Left 0
#define Right 1

/*******Macros*********/

#define idx(x) x
#define Tidx(var, x) (xDim * var + x)
