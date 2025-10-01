#ifndef CHOLESKYSOLVER_H_
#define CHOLESKYSOLVER_H_

#include <cmath>

void Cholesky_Decomposition(long double A[], int);

void Cholesky_BackSub(long double A[], int, int, long double B[]);

#endif // CHOLESKYSOLVER_H_
