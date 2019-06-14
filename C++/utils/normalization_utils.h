#ifndef GMMREG_UTILS_NORMALIZATION_UTILS_H_
#define GMMREG_UTILS_NORMALIZATION_UTILS_H_

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>

namespace gmmreg {

template <typename T>
void Normalize(vnl_matrix<T>& x, vnl_vector<T>& centroid, T& scale);

template <typename T>
void Denormalize(vnl_matrix<T>& x, const vnl_vector<T>& centroid,
                 const T scale);


}  // namespace gmmreg

#include "normalization_utils.cc"

#endif // GMMREG_UTILS_NORMALIZATION_UTILS_H_
