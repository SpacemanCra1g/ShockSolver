#ifndef DEFINITIONS_H_
#define DEFINITIONS_H_

/*********************************************/
/******* Space Reconstruction Methods ********/
/*********************************************/
#define WENO 1
#define FOG 2
#define GPR1 3
#define GPR2 4
#define MOOD531 5

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
#define SHOCKTUBE 2

/********************************/
/******* Indexing Macros ********/
/********************************/
#define DENS 0
#define MOMX 1
#define MOMY 2
#define MOMZ 3
#define ENER 4

#define DENSP 0
#define VELX 1
#define VELY 2
#define VELZ 3
#define PRES 4

#define LEFT 0
#define RIGHT 1

/*********************************/
/******* MISC Definitions ********/
/*********************************/
#define NDIMS 1
#define NumVar 5
#define NGC 3

#endif // DEFINITIONS_H_
