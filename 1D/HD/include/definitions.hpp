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
#define ENER 2

#define DENSP 0
#define VELX 1
#define PRES 2

#define LEFT 0
#define RIGHT 1

/*********************************/
/******* MISC Definitions ********/
/*********************************/
#define NDIMS 1
#define NumVar 3
#define NGC 3

#endif // DEFINITIONS_H_
