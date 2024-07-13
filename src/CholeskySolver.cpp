#include "../include/CholeskySolver.hpp"

void Cholesky_Decomposition(double *__restrict A, int n) {

  double sum1 = 0.0;
  double sum2 = 0.0;
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

void Cholesky_BackSub(double *__restrict A, int n, int m,
                      double *__restrict B) {

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
