#ifndef GMMREG_UTILS_GAUSS_TRANSFORM_H_
#define GMMREG_UTILS_GAUSS_TRANSFORM_H_

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>

namespace gmmreg {

template <typename T>
T GaussTransform(const vnl_matrix<T>& A, const vnl_matrix<T>& B, T scale);

template <typename T>
T GaussTransform(const vnl_matrix<T>& A, const vnl_matrix<T>& B, T scale,
                 vnl_matrix<T>& gradient);

}  // namespace gmmreg

#include "gauss_transform.cc"

#endif // GMMREG_UTILS_GAUSS_TRANSFORM_H_
