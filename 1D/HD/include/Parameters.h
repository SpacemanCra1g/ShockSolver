#include "definitions.hpp"

/****************
 * X Parameters *
 ****************/
#define NX 1000
#define X0 0.0
#define XN 9.0

/****************
 * T Parameters *
 ****************/
#define T0 0.0
#define TN 1.8

/*****************
 *Run Parameters *
 *****************/
#define SpaceMethod WENO
#define TestProblem SHUOSHER
#define BCS NEUMANN
#define CFL 0.8
#define EOS IDEAL
#define RK_Method 3
#define ell 6.0
#define MoodOrder 5
#define SlowStart false
#define GAMMA (1.4)
/****************
 *NN Parameters *
 ****************/
#define Use_NN false
#define UseDMP true
#define NegativeSlope 0.01
#define NN_Thresh .99

/******************************
 * Compiler defined constants *
 * Don't touch                *
 ******************************/
#define xDim (2 * NGC + NX)
#define XStart NGC
#define REdgeX xDim
#define XEnd (REdgeX - NGC)
#define dx ((XN - X0) / NX)

#define sigdel (ell * 1.4142135623730950488016887240L)

/*******Macros*********/
#define idx(x) x
#define Tidx(var, x) (xDim * var + x)
