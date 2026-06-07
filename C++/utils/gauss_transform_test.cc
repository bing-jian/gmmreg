#include <gtest/gtest.h>
#include <vnl/vnl_matrix.h>

#include "gauss_transform.h"

using namespace gmmreg;

// Two identical single points: exp(0) / 1 = 1.0
TEST(GaussTransform, IdenticalSinglePoints) {
    vnl_matrix<double> A(1, 2), B(1, 2);
    A(0, 0) = 1.0; A(0, 1) = 2.0;
    B(0, 0) = 1.0; B(0, 1) = 2.0;
    EXPECT_DOUBLE_EQ(GaussTransform(A, B, 1.0), 1.0);
}

// Points far apart with small scale: cross-term ≈ 0
TEST(GaussTransform, DistantPointsSmallScale) {
    vnl_matrix<double> A(1, 2), B(1, 2);
    A(0, 0) = 0.0;   A(0, 1) = 0.0;
    B(0, 0) = 100.0; B(0, 1) = 100.0;
    EXPECT_NEAR(GaussTransform(A, B, 0.01), 0.0, 1e-10);
}

// GaussTransform(A, B, s) == GaussTransform(B, A, s)
TEST(GaussTransform, Symmetry) {
    vnl_matrix<double> A(3, 2), B(4, 2);
    A(0, 0) = 0.0; A(0, 1) = 0.0;
    A(1, 0) = 1.0; A(1, 1) = 0.0;
    A(2, 0) = 0.0; A(2, 1) = 1.0;
    B(0, 0) = 0.5; B(0, 1) = 0.5;
    B(1, 0) = 1.5; B(1, 1) = 0.5;
    B(2, 0) = 0.5; B(2, 1) = 1.5;
    B(3, 0) = 1.5; B(3, 1) = 1.5;
    double scale = 0.5;
    EXPECT_NEAR(GaussTransform(A, B, scale), GaussTransform(B, A, scale), 1e-10);
}

// Larger scale → higher cross-term for the same point pair
TEST(GaussTransform, ScaleSensitivity) {
    vnl_matrix<double> A(1, 2), B(1, 2);
    A(0, 0) = 0.0; A(0, 1) = 0.0;
    B(0, 0) = 1.0; B(0, 1) = 0.0;
    EXPECT_LT(GaussTransform(A, B, 0.1), GaussTransform(A, B, 10.0));
}

// Result is always in (0, 1]
TEST(GaussTransform, ResultInRange) {
    vnl_matrix<double> A(3, 2), B(3, 2);
    A.fill(0); B.fill(0);
    A(0, 0) = 0.1; A(1, 0) = 0.5; A(2, 0) = 0.9;
    B(0, 1) = 0.2; B(1, 1) = 0.6; B(2, 1) = 0.8;
    double result = GaussTransform(A, B, 0.3);
    EXPECT_GT(result, 0.0);
    EXPECT_LE(result, 1.0);
}

// Verify gradient against finite differences
TEST(GaussTransform, GradientNumericalCheck) {
    vnl_matrix<double> A(3, 2), B(4, 2);
    A(0, 0) = 0.1; A(0, 1) = 0.2;
    A(1, 0) = 0.5; A(1, 1) = 0.3;
    A(2, 0) = 0.8; A(2, 1) = 0.7;
    B(0, 0) = 0.2; B(0, 1) = 0.4;
    B(1, 0) = 0.6; B(1, 1) = 0.1;
    B(2, 0) = 0.3; B(2, 1) = 0.9;
    B(3, 0) = 0.7; B(3, 1) = 0.5;
    double scale = 0.5;

    vnl_matrix<double> gradient(3, 2, 0.0);
    GaussTransform(A, B, scale, gradient);

    const double eps = 1e-5;
    for (int i = 0; i < 3; ++i) {
        for (int d = 0; d < 2; ++d) {
            vnl_matrix<double> A_plus = A, A_minus = A;
            A_plus(i, d) += eps;
            A_minus(i, d) -= eps;
            double fd = (GaussTransform(A_plus, B, scale) -
                         GaussTransform(A_minus, B, scale)) / (2.0 * eps);
            EXPECT_NEAR(gradient(i, d), fd, 1e-5)
                << "Gradient mismatch at (" << i << ", " << d << ")";
        }
    }
}
