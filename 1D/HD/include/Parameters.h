#include "definitions.hpp"

#define NDIMS 1
#define NGC 3

/****************
 * X Parameters *
 ****************/
#define NX 200
#define X0 0.0
#define XN 9.0

/****************
 * Y Parameters *
 ****************/
#define NY 1
#define Y0 0.0
#define YN 1.0
#define NumVar 3

/****************
 * T Parameters *
 ****************/
#define T0 0.0
#define TN 1.8

/*****************
 *Run Parameters *
 *****************/
#define SpaceMethod Weno

#define GAMMA 1.4
#define CFL 0.8
#define TestProblem ShuOsher
#define BoundaryCon Shu

#define RK_Method 3
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
#define Ener 2

#define Dens 0
#define VelX 1
#define Pres 2

#define Left 0
#define Right 1

/*******Macros*********/

#define idx(x) (x)
#define Tidx(var, x) (xDim * var + x)
