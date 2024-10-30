#include "../include/GP_Kernel.hpp"

long double GP_Kernel::int_egrand(long double x, double t) {

  return std::erf((t - x + .5) / sigdel) - std::erf((t - x - .5) / sigdel);
}

long double GP_Kernel::quad_cross(long double x, double t) {

  return (ell)*std::sqrt(.5 * M_PI) * int_egrand(x, t);
}

long double GP_Kernel::intg_predvec(long double x, int dir) {

  // These are flipped from the paper, and I think the github version is
  // wrong
  if (dir == LEFT) {
    return quad_cross(x, -0.5);
  } else {
    return quad_cross(x, 0.5);
  }
}

long double GP_Kernel::intg_kernel(long double x, long double y) {
  return quad_exact(x, y);
}

long double GP_Kernel::quad_exact(long double y, long double x) {
  long double yxp, yxn, yxm;
  yxp = (x - y + 1.0) / sigdel;
  yxn = (x - y) / sigdel;
  yxm = (x - y - 1.0) / sigdel;

  return std::sqrt(M_PI) * ell * ell *
         (yxp * std::erf(yxp) + yxm * std::erf(yxm) -
          2.L * (yxn * std::erf(yxn) +
                 (1.L / std::sqrt(M_PI)) * std::exp(-yxn * yxn)) +
          (1.L / std::sqrt(M_PI)) *
              (std::exp(-yxp * yxp) + std::exp(-yxm * yxm)));
}

void GP_Kernel::calculate_Preds1D(const int R) {
  const int StencilSize = 2 * R + 1;
  long double *stencil = new long double[StencilSize];
  long double *C = new long double[StencilSize * StencilSize];
  long double *TLeft = new long double[2 * StencilSize];
  long double *TRight = &TLeft[StencilSize];

  for (int i = 0; i < StencilSize; ++i) {
    stencil[i] = (long double)(i - R);
  }

  for (int i = 0; i < StencilSize; ++i) {
    for (int j = 0; j < StencilSize; ++j) {
      C[j + i * StencilSize] = intg_kernel(stencil[i], stencil[j]);
    }

    TRight[i] = intg_predvec(stencil[i], LEFT);
    TLeft[i] = intg_predvec(stencil[i], RIGHT);
  }

  Cholesky_Decomposition(C, StencilSize);

  Cholesky_BackSub(C, StencilSize, 2, TLeft);

  long double Sum = 0.0;

  for (int i = 0; i < StencilSize; ++i) {
    // Sum += std::fabs(TLeft[i]);
    Sum += TLeft[i];
  }
  for (int i = 0; i < StencilSize; ++i) {
    TLeft[i] /= Sum;
  }

  Sum = 0.0;
  for (int i = 0; i < StencilSize; ++i) {
    Sum += TRight[i];
    // Sum += std::fabs(TRight[i]);
  }
  for (int i = 0; i < StencilSize; ++i) {
    TRight[i] /= Sum;
  }

  // std::cout << std::endl;
  // for (int i = 0; i < StencilSize * 2; ++i) {
  //   std::cout << TLeft[i] << std::endl;
  //   if (i == StencilSize - 1) {
  //     std::cout << std::endl;
  //   }
  // }

  if (R == 1) {
    // R1Left[0] = (double)TLeft[1];
    // R1Left[1] = (double)TLeft[0];
    // R1Left[2] = (double)TLeft[2];
    // R1Right[0] = (double)TRight[1];
    // R1Right[1] = (double)TRight[0];
    // R1Right[2] = (double)TRight[2];
    for (int i = 0; i < 3; ++i) {
      R1Left[i] = (double)TLeft[i];
      R1Right[i] = (double)TRight[i];
    }
  } else if (R == 2) {
    for (int i = 0; i < 5; ++i) {
      R2Left[i] = (double)TLeft[i];
      R2Right[i] = (double)TRight[i];
    }
  }

  delete[] stencil;
  delete[] C;
  delete[] TLeft;
}
