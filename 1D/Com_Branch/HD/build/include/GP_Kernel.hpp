#ifndef GP_KERNEL_H_
#define GP_KERNEL_H_

#include "CholeskySolver.hpp"
#include "SourceParameters.h"
#include <cblas.h>
#include <cmath>
#include <iostream>

class GP_Kernel {

private:
  long double int_egrand(long double x, double t);

  long double quad_cross(long double x, double t);

  long double intg_predvec(long double x, int dir);

  long double intg_kernel(long double x, long double y);

  long double quad_exact(long double y, long double x);

public:
  double R1Left[3];
  double R1Right[3];

  double R2Left[5];
  double R2Right[5];

  void calculate_Preds1D(const int R);
};

#endif // GP_KERNEL_H_
