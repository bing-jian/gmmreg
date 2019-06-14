#include "normalization_utils.h"

namespace gmmreg {

template <typename T>
void Normalize(vnl_matrix<T>& x, vnl_vector<T>& centroid, T& scale) {
  int n = x.rows();
  if (n == 0) return;
  int d = x.cols();
  centroid.set_size(d);

  vnl_vector<T> col;
#pragma omp for
  for (int i = 0; i < d; ++i) {
    col = x.get_column(i);
    centroid(i) = col.mean();
  }
#pragma omp for
  for (int i = 0; i < n; ++i) {
    x.set_row(i, x.get_row(i) - centroid);
  }
  scale = x.frobenius_norm() / sqrt(T(n));
  x = x / scale;
}

template <typename T>
void Denormalize(vnl_matrix<T>& x, const vnl_vector<T>& centroid,
                 const T scale) {
  int n = x.rows();
  if (n == 0) return;
  int d = x.cols();
#pragma omp for
  for (int i = 0; i < n; ++i) {
    x.set_row(i, x.get_row(i) * scale + centroid);
  }
}

}  // namespace gmmreg
