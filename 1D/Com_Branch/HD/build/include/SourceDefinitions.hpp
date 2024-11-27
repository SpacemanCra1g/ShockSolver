#ifndef DEFINITIONS_H_
#define DEFINITIONS_H_

/*********************************************/
/******* Space Reconstruction Methods ********/
/*********************************************/
#define WENO 1
#define FOG 2
#define GPR1 3
#define GPR2 4
#define MOOD 5
#define PLM 6

/***********************************/
/******* Equations of State ********/
/***********************************/
#define IDEAL 1

/************************************/
/******* Boundary Conditions ********/
/************************************/
#define NEUMANN 1
#define SHOCKTUBEBC 2

/**********************************/
/******* Inital Conditions ********/
/**********************************/
#define SHUOSHER 1
#define SRSHOCKTUBE 2

/*********************************/
/******* Riemann Solvers ********/
/*********************************/
#define HLL 1
#define HLLC 2

/*******************************/
/******* Slope Limiters ********/
/*******************************/
#define VANLEER 1
#define MC 2
#define MINMOD 3

/********************************/
/******* Indexing Macros ********/
/********************************/
#define LEFT 0
#define RIGHT 1

#endif // DEFINITIONS_H_
