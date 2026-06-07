#include <vector>

#include <gtest/gtest.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>

#include "rotation_utils.h"

using namespace gmmreg;

static vnl_matrix<double> mat33() { return vnl_matrix<double>(3, 3, 0.0); }

static vnl_vector<double> quat(double x, double y, double z, double r) {
    vnl_vector<double> q(4);
    q[0]=x; q[1]=y; q[2]=z; q[3]=r;
    return q;
}

// ── Quaternion2Rotation (simple overload) ────────────────────────────────────

// q=[0,0,0,1] (identity rotation) → R = I
TEST(Quaternion2Rotation, IdentityQuaternion) {
    vnl_matrix<double> R = mat33();
    Quaternion2Rotation(quat(0, 0, 0, 1), R);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            EXPECT_NEAR(R(i, j), (i == j ? 1.0 : 0.0), 1e-10);
}

// q=[1,0,0,0] (180° around x) → diag(1, -1, -1)
TEST(Quaternion2Rotation, HalfTurnAroundX) {
    vnl_matrix<double> R = mat33();
    Quaternion2Rotation(quat(1, 0, 0, 0), R);
    EXPECT_NEAR(R(0, 0),  1.0, 1e-10);
    EXPECT_NEAR(R(1, 1), -1.0, 1e-10);
    EXPECT_NEAR(R(2, 2), -1.0, 1e-10);
    EXPECT_NEAR(R(0, 1),  0.0, 1e-10);
    EXPECT_NEAR(R(1, 0),  0.0, 1e-10);
}

// R must be orthogonal: R * R^T = I
TEST(Quaternion2Rotation, Orthogonality) {
    vnl_matrix<double> R = mat33();
    Quaternion2Rotation(quat(0.5, 0.5, 0.5, 0.5), R);
    vnl_matrix<double> RRt = R * R.transpose();
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            EXPECT_NEAR(RRt(i, j), (i == j ? 1.0 : 0.0), 1e-10);
}

// Scaling q by a constant must not change R (normalised by ||q||^2)
TEST(Quaternion2Rotation, ScaleInvariance) {
    vnl_matrix<double> R1 = mat33(), R2 = mat33();
    auto q = quat(1, 2, 3, 4);
    Quaternion2Rotation(q,       R1);
    Quaternion2Rotation(q * 5.0, R2);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            EXPECT_NEAR(R1(i, j), R2(i, j), 1e-10);
}

// ── Quaternion2Rotation with gradients ──────────────────────────────────────

// g1..g4 must match central finite differences of R w.r.t. q[0..3].
// The gradient formulas in rotation_utils.cc are only correct for unit
// quaternions (ss=1 cancels the ss vs ss² difference in off-diagonal terms).
TEST(Quaternion2Rotation, GradientNumericalCheck) {
    auto q0 = quat(1, 2, 3, 4);
    q0 = q0 / q0.two_norm();  // normalise: gradient code requires ||q||=1
    vnl_matrix<double> R = mat33();
    vnl_matrix<double> g1 = mat33(), g2 = mat33(), g3 = mat33(), g4 = mat33();
    Quaternion2Rotation(q0, R, g1, g2, g3, g4);

    const double eps = 1e-6;
    std::vector<vnl_matrix<double>*> grads = {&g1, &g2, &g3, &g4};
    for (int k = 0; k < 4; ++k) {
        auto qp = q0, qm = q0;
        qp[k] += eps;
        qm[k] -= eps;
        vnl_matrix<double> Rp = mat33(), Rm = mat33();
        Quaternion2Rotation(qp, Rp);
        Quaternion2Rotation(qm, Rm);
        vnl_matrix<double> fd = (Rp - Rm) / (2.0 * eps);
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                EXPECT_NEAR((*grads[k])(i, j), fd(i, j), 1e-5)
                    << "g" << k+1 << "(" << i << "," << j << ") mismatch";
    }
}
