#include "../include/CholeskySolver.hpp"
#include <iostream>

void Cholesky_Decomposition(long double A[], int n) {

  long double sum1 = 0.0;
  long double sum2 = 0.0;
  // Decomposing a matrix into Lower Triangular
  for (int i = 0; i < n; ++i) {
    sum1 = 0.0;
    for (int k = 0; k < i; k++) {
      sum1 += A[k * n + i] * A[k * n + i];
    }
    A[i * n + i] = std::sqrt(A[i * n + i] - sum1);

    for (int j = i + 1; j < n; ++j) {
      sum2 = 0.0;
      for (int k = 0; k < i; ++k) {
        sum2 += A[k * n + j] * A[k * n + i];
      }
      A[i * n + j] = (A[i * n + j] - sum2) / A[i * n + i];
    }
  }
}

// void Cholesky_Decomposition(long double A[], int n) {
//   long double sum1 = 0.0;
//   long double sum2 = 0.0;

//   // Decomposing a matrix into Lower Triangular
//   for (int i = 0; i < n; ++i) {
//     sum1 = 0.0;

//     // Compute diagonal elements
//     for (int k = 0; k < i; ++k) {
//       sum1 += A[i * n + k] *
//               A[i * n + k]; // Note: A[i * n + k] instead of A[k * n + i]
//     }

//     long double diag_val = A[i * n + i] - sum1;

//     if (diag_val <= 0.0) {
//       throw std::runtime_error("Matrix is not positive definite.");
//     }

//     A[i * n + i] = std::sqrt(diag_val);

//     // Compute off-diagonal elements
//     for (int j = i + 1; j < n; ++j) {
//       sum2 = 0.0;
//       for (int k = 0; k < i; ++k) {
//         sum2 += A[j * n + k] *
//                 A[i * n + k]; // Note: A[j * n + k] instead of A[k * n + i]
//       }
//       A[j * n + i] = (A[j * n + i] - sum2) / A[i * n + i];
//     }
//   }
// }

void Cholesky_BackSub(long double A[], int n, int m, long double B[]) {

  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) {
      for (int k = j - 1; k >= 0; --k) {
        B[i * n + j] -= B[i * n + k] * A[k * n + j];
      }
      B[i * n + j] /= A[j * n + j];
    }

    for (int j = n - 1; j >= 0; --j) {
      for (int k = j + 1; k < n; ++k) {
        B[i * n + j] -= B[i * n + k] * A[k + n * j];
      }
      B[i * n + j] /= A[j * n + j];
    }
  }
}
