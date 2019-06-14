#include "em_utils.h"

namespace gmmreg {

template <typename T>
void ComputeP(const vnl_matrix<T>& x, const vnl_matrix<T>& y, vnl_matrix<T>& P,
              T& E, T sigma, int outliers) {
  T k;
  k = -2 * sigma * sigma;

  vnl_vector<T> column_sum;
  int m = x.rows();
  int s = y.rows();
  int d = x.cols();
  column_sum.set_size(s);
  column_sum.fill(0);
  T outlier_term = outliers * pow((2 * sigma * sigma * 3.1415926), 0.5 * d);
#pragma omp for
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < s; ++j) {
      T r = 0;
      for (int t = 0; t < d; ++t) {
        r += (x(i, t) - y(j, t)) * (x(i, t) - y(j, t));
      }
      P(i, j) = exp(r / k);
      column_sum[j] += P(i, j);
    }
  }

  if (outliers != 0) {
#pragma omp for
    for (int i = 0; i < s; ++i) column_sum[i] += outlier_term;
  }
  if (column_sum.min_value() > (1e-12)) {
    E = 0;
#pragma omp for
    for (int i = 0; i < s; ++i) {
      for (int j = 0; j < m; ++j) {
        P(j, i) = P(j, i) / column_sum[i];
      }
      E -= log(column_sum[i]);
    }
  } else {
    P.empty();
  }
}

}  // namespace gmmreg
