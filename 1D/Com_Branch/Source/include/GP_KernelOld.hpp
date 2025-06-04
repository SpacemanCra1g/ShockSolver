/*
#ifndef GP_KERNEL_H_
#define GP_KERNEL_H_

#define YDIR 1
#define XDIR 2

#define Rg1 1
#define Rg2 2

#define Tg1 3
#define Tg2 4

#define Bg1 5
#define Bg2 6

#define Lg1 7
#define Lg2 8
#define sth 1.0 / (2.0 * std::sqrt(3.0))

#include "CholeskySolver.hpp"
#include "Parameters.h"
#include <cblas.h>
#include <cmath>
#include <iostream>

class GP_Kernel {

private:
  double R1_Vec1D[2 * 3];
  double R1_VecX1D[2 * 3] = {1.5, .5, -.5, .5, -.5, -1.5};
  double R1_C1D[9];
  double R1_DeltX1D[9] = {0.0, 1.0, 2.0, -1.0, 0.0, 1.0, -2.0, -1.0, 0.0};

  // Underlying structure of R1 prediction vectors
  double R1_Vec[8 * 5]; // 8 vectors length 5
                        // This is for 2 quad points on each face

  double R1_VecX[40] = { // one length 5 vec for each quad point
                         // reconstruction in the X direction
      .5,         1.5,     .5,         -.5,        .5,   .5,         1.5,
      .5,         -.5,     .5,         -.5,        .5,   -.5,        -1.5,
      -.5,        -.5,     .5,         -.5,        -1.5, -.5,        sth,
      1 + sth,    sth,     -1.0 + sth, sth,        -sth, 1 - sth,    -sth,
      -1.0 - sth, -sth,    sth,        1 + sth,    sth,  -1.0 + sth, sth,
      -sth,       1 - sth, -sth,       -1.0 - sth, -sth};

  double R1_VecY[40] = { // one length 5 vec for each quad point
                         // reconstruction in the Y direction
      sth, sth, -1 + sth, sth, 1 + sth, -sth, -sth, -1 - sth, -sth, 1 - sth,
      sth, sth, -1 + sth, sth, 1 + sth, -sth, -sth, -1 - sth, -sth, 1 - sth,
      .5,  .5,  -.5,      .5,  1.5,     .5,   .5,   -.5,      .5,   1.5,
      -.5, -.5, -1.5,     -.5, .5,      -.5,  -.5,  -1.5,     -.5,  .5};

  double R1_C[25]; // Radius 1 2D kernal

  // index by (m-1)*5 + (n-1)
  // These are the offset ammounts for the 2D 2qp
  double R1_DeltY[25] = {0.0,  0.0,  1.0,  0.0, -1.0, 0.0,  0.0, 1.0, 0.0,
                         -1.0, -1.0, -1.0, 0.0, -1.0, -2.0, 0.0, 0.0, 1.0,
                         0.0,  -1.0, 1.0,  1.0, 2.0,  1.0,  0.0};

  double R1_DeltX[25] = {0.0, -1.0, 0.0,  1.0,  0.0, 1.0, 0.0,  1.0,  2.0,
                         1.0, 0.0,  -1.0, 0.0,  1.0, 0.0, -1.0, -2.0, -1.0,
                         0.0, -1.0, 0.0,  -1.0, 0.0, 1.0, 0.0};

public:
  // These are the public access points for the
  // Gaussian quadrature prediction vectors
  double *R1_g1R = &R1_Vec[0];
  double *R1_g2R = &R1_Vec[5];
  double *R1_g1L = &R1_Vec[10];
  double *R1_g2L = &R1_Vec[15];
  double *R1_g1U = &R1_Vec[20];
  double *R1_g2U = &R1_Vec[25];
  double *R1_g1D = &R1_Vec[30];
  double *R1_g2D = &R1_Vec[35];

  double *R1_g1R1D = R1_Vec1D;
  double *R1_g1L1D = R1_Vec1D + 3;

  void GP_Kernel_init() {

    const double valX = std::sqrt(2.0) * ell / dx;
    const double valY = std::sqrt(2.0) * ell / dy;

    double ax, bx, cx, ay, by, cy;

    for (int i = 0; i < 25; ++i) {
      ax = (R1_DeltX[i] + 1) / valX;
      bx = (R1_DeltX[i] - 1) / valX;
      cx = (R1_DeltX[i]) / valX;

      ay = (R1_DeltY[i] + 1) / valY;
      by = (R1_DeltY[i] - 1) / valY;
      cy = (R1_DeltY[i]) / valY;

      R1_C[i] =
          (M_PI * valX * valY / 4.0) *
          (ax * std::erf(ax) + bx * std::erf(bx) +
           (1 / std::sqrt(M_PI)) * (std::exp(-ax * ax) + std::exp(-bx * bx)) -
           2 * (cx * std::erf(cx) +
                (1 / std::sqrt(M_PI)) * std::exp(-cx * cx))) *
          (ay * std::erf(ay) + by * std::erf(by) +
           (1 / std::sqrt(M_PI)) * (std::exp(-ay * ay) + std::exp(-by * by)) -
           2 * (cy * std::erf(cy) +
                (1 / std::sqrt(M_PI)) * std::exp(-cy * cy)));
    }

    for (int i = 0; i < 40; ++i) {
      ax = (R1_VecX[i] + .5) / valX;
      bx = (R1_VecX[i] - .5) / valX;

      ay = (R1_VecY[i] + .5) / valY;
      by = (R1_VecY[i] - .5) / valY;

      R1_Vec[i] = M_PI * ell * ell / (2.0 * dx * dy) *
                  (std::erf(ax) - std::erf(bx)) * (std::erf(ay) - std::erf(by));
    }

    Cholesky_Decomposition(R1_C, 5);
    Cholesky_BackSub(R1_C, 5, 8, R1_Vec);

    double Sum1;
    for (int k = 0; k < 8; ++k) {
      Sum1 = 0.0;
      for (int i = 0; i < 5; ++i) {
        Sum1 += R1_Vec[5 * k + i];
      }

      for (int i = 0; i < 5; ++i) {
        R1_Vec[5 * k + i] /= Sum1;
      }
    }

    cblas_dscal(40, 0.5, R1_Vec, 1);
  }

  void TVector(int ndims, int NQP, int qp, int dir, int vecSize,
               double *inputVec) {
    switch (ndims) {
    case 1:
      switch (dir) {
      case Right:
        switch (vecSize) {}
      }
    }
  }
};

function intg_kernel(x, y) result(f)
    implicit none
    real, intent(IN) :: x, y
    real :: f
    f = quad_exact(x, y)
    return
end function intg_kernel


function quad_exact(x1,x2) result(Integ)

    !exact quadrature, only good for SE kernel
    real, intent(IN) :: x1, x2

    real :: Integ, yxp, yxn, yxm, sigdel
    sigdel = sim_sigdel*SQRT(2.)

    yxp = (x1 - x2 + 1.)/sigdel
    yxn = (x1      -x2)/sigdel
    yxm = (x1 - x2 -1.)/sigdel


    Integ = 0.5*SQRT(PI)*(sigdel)**2 *( yxp*ERF(yxp) + yxm*ERF(yxm) &
         - 2.*( yxn*ERF(yxn) + 1./SQRT(PI) *EXP(-yxn**2) ) &
         + 1./SQRT(PI) * ( EXP(-yxp**2) + exp(-yxm**2) ) )
    return
end function quad_exact

////
////
    do i = 1, N
     stencil(i) =  REAL(i - R - 1)
  end do

  !first thing is to calculate the covariance matrix according to eq. 15
  !since C is symmetric only bother with one side of the diaganol
  do i = 1,N
     do j = 1,N
        C(i,j) = intg_kernel(stencil(i), stencil(j))
     end do
    T(1:2, i) = intg_predvec(stencil(i))
  end do
  !now we need to solve the linear eqns for v & Z (see eqs 30-32)

#endif // GP_KERNEL_H_
*/
