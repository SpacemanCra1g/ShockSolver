#ifndef CHOLESKYSOLVER_H_
#define CHOLESKYSOLVER_H_

#include <cmath>

void Cholesky_Decomposition(double *__restrict, int);

void Cholesky_BackSub(double *__restrict, int, int, double *__restrict);

#endif // CHOLESKYSOLVER_H_
