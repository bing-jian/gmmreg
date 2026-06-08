#include "transform_utils.h"

#include <algorithm>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace gmmreg {

vnl_matrix<double> RelativeTransform(const vnl_matrix<double>& T_src,
                                     const vnl_matrix<double>& T_tgt) {
  const vnl_matrix<double> R_src = T_src.extract(3, 3, 0, 0);
  const vnl_matrix<double> R_tgt = T_tgt.extract(3, 3, 0, 0);
  const vnl_vector<double> t_src(T_src.get_column(3).extract(3));
  const vnl_vector<double> t_tgt(T_tgt.get_column(3).extract(3));

  const vnl_matrix<double> R_rel = R_tgt.transpose() * R_src;
  const vnl_vector<double> t_rel = R_tgt.transpose() * (t_src - t_tgt);

  vnl_matrix<double> T(4, 4, 0.0);
  T.update(R_rel, 0, 0);
  for (int i = 0; i < 3; ++i) T(i, 3) = t_rel[i];
  T(3, 3) = 1.0;
  return T;
}

double RotationErrorDeg(const vnl_matrix<double>& T_est,
                        const vnl_matrix<double>& T_ref) {
  const vnl_matrix<double> R_est = T_est.extract(3, 3, 0, 0);
  const vnl_matrix<double> R_ref = T_ref.extract(3, 3, 0, 0);
  const vnl_matrix<double> R     = R_ref.transpose() * R_est;
  double cos_theta = (R(0, 0) + R(1, 1) + R(2, 2) - 1.0) * 0.5;
  cos_theta = std::max(-1.0, std::min(1.0, cos_theta));
  return std::acos(cos_theta) * 180.0 / M_PI;
}

double RotationMagnitudeDeg(const vnl_matrix<double>& T) {
  const vnl_matrix<double> R = T.extract(3, 3, 0, 0);
  double cos_theta = (R(0, 0) + R(1, 1) + R(2, 2) - 1.0) * 0.5;
  cos_theta = std::max(-1.0, std::min(1.0, cos_theta));
  return std::acos(cos_theta) * 180.0 / M_PI;
}

}  // namespace gmmreg
