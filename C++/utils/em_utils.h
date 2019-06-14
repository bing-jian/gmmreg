#ifndef GMMREG_UTILS_EM_UTILS_H_
#define GMMREG_UTILS_EM_UTILS_H_

#include <vnl/vnl_matrix.h>

namespace gmmreg {

template <typename T>
void ComputeP(const vnl_matrix<T>& x, const vnl_matrix<T>& y, vnl_matrix<T>& P,
              T& E, T sigma, int outliers);

}  // namespace gmmreg

#include "em_utils.cc"

#endif // GMMREG_UTILS_EM_UTILS_H_
