#define NDIMS 1
#define NGC 3

/***************
 * X Parameters *
 ****************/
#define NX 10
#define X0 0
#define XN 9

/***************
 * Y Parameters *
 ****************/
#if NDIMS == 1
#define NY 1
#define Y0 0
#define YN 1
#define NumVar 3

#else
#define NY 1
#define Y0 0
#define YN 9
#define NumVar 4

#endif

/***************
 * T Parameters *
 ****************/
#define T0 0
#define TN 2

/****************
 *Run Parameters *
 *****************/

#define GAMMA 1.4
#define CFL 0.8
#define TestProblem ShuOsher
#define BC ShuOsher

#define RK_Method 1
#define ell 12
#define MoodOrder 5
#define SlowStart True

#define nqp 1

/***************
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
#define xDim 2 * NGC + NX
#define XStart NGC
#define REdgeX xDim
#define XEnd REdgeX - NGC
#define dx (XN - X0) / NX

#if NDIMS == 1
#define yDim 1
#define YStart 0
#define REdgeY 1
#define YEnd 1
#define dy dx

#else
#define yDim 2 * NGC + NY
#define YStart NGC
#define REdgeY yDim
#define YEnd REdgeY - NGC
#define dy (YN - Y0) / NY
#endif

#define Dens 0
#define MomX 1
#define Ener 2
#define MomY 3

#define Left 0
#define Right 1
#define Bottom 2
#define Top 3

/*******Macros*********/
#if NDIMS == 2
#define idx(x, y) x *yDim + y
#else
#define idx(x, y) x
#endif
