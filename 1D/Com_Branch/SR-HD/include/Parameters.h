#ifndef PARAM_H_
#define PARAM_H_

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

#define PL 1000.0
#define PR 0.01

#define XVELL 0.0
#define YVELL 0.99
#define ZVELL 0.0

#define XVELR 0.0
#define YVELR 0.0
#define ZVELR 0.0

/****************
 * T Parameters *
 ****************/
#define T0 0.0
#define TN 0.4

/*****************
 *Run Parameters *
 *****************/
#define EvolveChars false
#define SpaceMethod MOOD
#define TestProblem SHOCKTUBE
#define BCS SHOCKTUBEBC
#define RIEMANN HLLC
#define CFL 0.8
#define EOS IDEAL
#define RK_Method 3 /*set to -1 for CharTracing, sets automatically for PLM*/
#define ell 6.0
#define MoodOrder 3
#define SlowStart true
#define GAMMA (5.0 / 3.0)
#define LIMITSLOPE MINMOD
/****************
 *NN Parameters *
 ****************/
#define Use_NN false
#define UseDMP true
#define NegativeSlope 0.01
#define NN_Thresh .99

#endif
