#include "rbf_utils.h"

#include "macros.h"

namespace gmmreg {

template <typename T>
void GaussianAffinityMatrix(const T* A, const T* B, int m, int n, int dim,
                            T scale, T* dist) {
  scale = -2.0 * SQR(scale);
  int k = 0;
#pragma omp for
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) {
      T dist_ij = 0;
      for (int d = 0; d < dim; ++d) {
        dist_ij += SQR(A[i * dim + d] - B[j * dim + d]);
      }
      dist[k] = exp(dist_ij / scale);
      ++k;
    }
  }
}

// TODO: add one more version when the model is same as ctrl_pts
// reference:  Landmark-based Image Analysis, Karl Rohr, p195
template <typename T>
void ComputeTPSKernel(const vnl_matrix<T>& model, const vnl_matrix<T>& ctrl_pts,
                      vnl_matrix<T>& U, vnl_matrix<T>& K) {
  int m = model.rows();
  int n = ctrl_pts.rows();
  int d = ctrl_pts.cols();
  // asssert(model.cols()==d==(2|3));
  K.set_size(n, n);
  K.fill(0);
  U.set_size(m, n);
  U.fill(0);
  T eps = 1e-006;

#pragma omp for
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) {
      vnl_vector<T> v_ij = model.get_row(i) - ctrl_pts.get_row(j);
      if (d == 2) {
        T r = v_ij.squared_magnitude();
        if (r > eps) {
          U(i, j) = r * log(r) / 2;
        }
      } else if (d == 3) {
        T r = v_ij.two_norm();
        U(i, j) = -r;
      }
    }
  }

#pragma omp for
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      vnl_vector<T> v_ij = ctrl_pts.get_row(i) - ctrl_pts.get_row(j);
      if (d == 2) {
        T r = v_ij.squared_magnitude();
        if (r > eps) {
          K(i, j) = r * log(r) / 2;
        }
      } else if (d == 3) {
        T r = v_ij.two_norm();
        K(i, j) = -r;
      }
    }
  }

#pragma omp for
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < i; ++j) {
      K(i, j) = K(j, i);
    }
  }
}

/*
   Matlab code in cpd_G.m:
   k=-2*beta^2;
   [n, d]=size(x);  [m, d]=size(y);

   G=repmat(x,[1 1 m])-permute(repmat(y,[1 1 n]),[3 2 1]);
   G=squeeze(sum(G.^2,2));
   G=G/k;
   G=exp(G);
   */
template <typename T>
void ComputeGaussianKernel(const vnl_matrix<T>& model,
                           const vnl_matrix<T>& ctrl_pts, vnl_matrix<T>& G,
                           vnl_matrix<T>& K, T lambda) {
  int m, n, d;
  m = model.rows();
  n = ctrl_pts.rows();
  d = ctrl_pts.cols();
  // asssert(model.cols()==d);
  // assert(lambda>0);

  G.set_size(m, n);
  GaussianAffinityMatrix(model.data_block(), ctrl_pts.data_block(), m, n, d,
                         lambda, G.data_block());

  if (model == ctrl_pts) {
    K = G;
  } else {
    K.set_size(n, n);
    GaussianAffinityMatrix(ctrl_pts.data_block(), ctrl_pts.data_block(), n, n,
                           d, lambda, K.data_block());
  }
}

}  // namespace gmmreg
