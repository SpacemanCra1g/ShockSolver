#ifndef PARAM_H_
#define PARAM_H_

#include "definitions.hpp"

/****************
 * X Parameters *
 ****************/
#define NX 400
#define X0 0.0
#define XN 1.0

/****************
 * T Parameters *
 ****************/
#define T0 0.0
#define TN .25

/*****************
 *Run Parameters *
 *****************/
#define EvolveChars false
#define SpaceMethod MOOD
#define TestProblem SHOCKTUBE
#define BCS SHOCKTUBEBC
#define RIEMANN HLL
#define CFL 0.8
#define EOS IDEAL
#define RK_Method 3 /*set to -1 for CharTracing, sets automatically for PLM*/
#define ell 6.0
#define MoodOrder 5
#define SlowStart true
#define GAMMA (1.4)
#define LIMITSLOPE MINMOD
/****************
 *NN Parameters *
 ****************/
#define Use_NN false
#define UseDMP true
#define NegativeSlope 0.01
#define NN_Thresh .99

#endif
