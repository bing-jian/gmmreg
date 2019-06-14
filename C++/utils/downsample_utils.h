#ifndef GMMREG_UTILS_DOWNSAMPLE_UTILS_H_
#define GMMREG_UTILS_DOWNSAMPLE_UTILS_H_

#include <vnl/vnl_matrix.h>
#include "../fgt_utils.h"

namespace gmmreg {

template <typename T>
vnl_matrix<T> CreateSparseNodes(const NanoflannTree<T>& tree, T dist_radius);

}  // namespace gmmreg

#include "downsample_utils.cc"

#endif // GMMREG_UTILS_DOWNSAMPLE_UTILS_H_
