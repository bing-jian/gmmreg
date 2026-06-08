#ifndef GMMREG_UTILS_TRANSFORM_UTILS_H_
#define GMMREG_UTILS_TRANSFORM_UTILS_H_

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>

namespace gmmreg {

// Relative rigid transform T_tgt^{-1} * T_src.
// Both inputs are (d+1)×(d+1) homogeneous matrices with orthonormal R block.
// Uses R^{-1} = R^T (no general matrix inverse needed).
vnl_matrix<double> RelativeTransform(const vnl_matrix<double>& T_src,
                                     const vnl_matrix<double>& T_tgt);

// Geodesic rotation error in degrees between two 4×4 rigid transforms.
// Returns arccos((trace(R_ref^T * R_est) - 1) / 2) in degrees.
double RotationErrorDeg(const vnl_matrix<double>& T_est,
                        const vnl_matrix<double>& T_ref);

// Rotation magnitude (geodesic distance from identity) of a 4×4 rigid
// transform in degrees.
double RotationMagnitudeDeg(const vnl_matrix<double>& T);

}  // namespace gmmreg

#endif  // GMMREG_UTILS_TRANSFORM_UTILS_H_
