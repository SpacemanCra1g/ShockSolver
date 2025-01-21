#ifndef SRC_PARAM_H_
#define SRC_PARAM_H_

/* #include "definitions.hpp" */
#include "Parameters.h" // This just makes the LSP Work
#include "SourceDefinitions.hpp"
#include <cmath>

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

#endif
