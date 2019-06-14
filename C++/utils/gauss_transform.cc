#include "gauss_transform.h"

#include "macros.h"

namespace gmmreg {

/*
 *  Note: The input point set containing 'n' points in 'd'-dimensional
 *  space should be arranged in the memory such that the j-th coordinate of i-th
 * point
 *  is at location  (i*d + j), i.e., the input matrix in MATLAB should
 *  be of size (d, n) since MATLAB is column-major while the input array
 *  in NumPy should be of size (n, d) since NumPy is row-major.
 */
template <typename T>
T GaussTransform(const T* A, const T* B, int m, int n, int dim, T scale) {
  T cross_term = 0;
  scale = SQR(scale);
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) {
      T dist_ij = 0;
      for (int d = 0; d < dim; ++d) {
        dist_ij += SQR(A[i * dim + d] - B[j * dim + d]);
      }
      T cost_ij = exp(-1.0 * dist_ij / scale);
      cross_term += cost_ij;
    }
    /* printf("cross_term = %.3f\n", cross_term);  */
  }
  return cross_term / (m * n);
}

template <typename T>
T GaussTransform(const T* A, const T* B, int m, int n, int dim, T scale,
                 T* grad) {
  T cross_term = 0;

#pragma omp parallel for
  for (int i = 0; i < m * dim; ++i) {
    grad[i] = 0;
  }

  scale = SQR(scale);

#pragma omp parallel for reduction(+ : cross_term)
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) {
      T dist_ij = 0;
      for (int d = 0; d < dim; ++d) {
        dist_ij += SQR(A[i * dim + d] - B[j * dim + d]);
      }
      T cost_ij = exp(-1.0 * dist_ij / scale);
      for (int d = 0; d < dim; ++d) {
        grad[i * dim + d] -= cost_ij * 2.0 * (A[i * dim + d] - B[j * dim + d]);
      }
      cross_term += cost_ij;
    }
  }

  scale *= m * n;
#pragma omp parallel for
  for (int i = 0; i < m * dim; ++i) {
    grad[i] /= scale;
  }
  return cross_term / (m * n);
}

template <typename T>
T GaussTransform(const vnl_matrix<T>& A, const vnl_matrix<T>& B, T scale) {
  // assert A.cols() == B.cols()
  return GaussTransform(A.data_block(), B.data_block(), A.rows(), B.rows(),
                        A.cols(), scale);
}

template <typename T>
T GaussTransform(const vnl_matrix<T>& A, const vnl_matrix<T>& B, T scale,
                 vnl_matrix<T>& gradient) {
  // assert A.cols() == B.cols()
  return GaussTransform(A.data_block(), B.data_block(), A.rows(), B.rows(),
                        A.cols(), scale, gradient.data_block());
}

}  // namespace gmmreg
