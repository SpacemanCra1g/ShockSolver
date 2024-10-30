#include "definitions.hpp"

/****************
 * X Parameters *
 ****************/
#define NX 400
#define X0 0.0
#define XN 1.0

/*****************
 * ShockTube ICs *
 *****************/
#define RHOL 1.0
#define RHOR 1.0

#define PL 1000
#define PR 0.01

#define XVELL 0.0
#define YVELL 0.0
#define ZVELL 0.0

#define XVELR 0.0
#define YVELR 0.99
#define ZVELR 0.0

/****************
 * T Parameters *
 ****************/
#define T0 0.0
#define TN 0.4

/*****************
 *Run Parameters *
 *****************/
#define SpaceMethod WENO
#define TestProblem SHOCKTUBE
#define BCS SHOCKTUBEBC
#define CFL 0.8
#define EOS IDEAL
#define RK_Method 3
#define ell 6.0
#define MoodOrder 5
#define SlowStart true
#define GAMMA (5.0 / 3.0)
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
