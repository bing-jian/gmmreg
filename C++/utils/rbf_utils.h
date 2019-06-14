#ifndef GMMREG_UTILS_RBF_UTILS_H_
#define GMMREG_UTILS_RBF_UTILS_H_

#include <vnl/vnl_matrix.h>

namespace gmmreg {

template <typename T>
void GaussianAffinityMatrix(const T* A, const T* B, int m, int n, int dim,
                            T scale, T* dist);

template <typename T>
void ComputeTPSKernel(const vnl_matrix<T>& model, const vnl_matrix<T>& ctrl_pts,
                      vnl_matrix<T>& U, vnl_matrix<T>& K);

template <typename T>
void ComputeGaussianKernel(const vnl_matrix<T>& model,
                           const vnl_matrix<T>& ctrl_pts, vnl_matrix<T>& G,
                           vnl_matrix<T>& K, T beta);

}  // namespace gmmreg

#include "rbf_utils.cc"

#endif // GMMREG_UTILS_RBF_UTILS_H_
