#include <cmath>

#include <gtest/gtest.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>

#include "normalization_utils.h"

using namespace gmmreg;

// ── Normalize ────────────────────────────────────────────────────────────────

// Each column of x must have zero mean after normalization
TEST(Normalize, ZeroCentroidAfterNormalization) {
    vnl_matrix<double> x(4, 2);
    x(0, 0)=1; x(0, 1)=2;
    x(1, 0)=3; x(1, 1)=4;
    x(2, 0)=5; x(2, 1)=6;
    x(3, 0)=7; x(3, 1)=8;
    vnl_vector<double> centroid;
    double scale;
    Normalize(x, centroid, scale);
    for (int d = 0; d < 2; ++d) {
        double mean = 0;
        for (int i = 0; i < 4; ++i) mean += x(i, d);
        mean /= 4;
        EXPECT_NEAR(mean, 0.0, 1e-10) << "column " << d << " not zero-centred";
    }
}

// ||x||_F / sqrt(n) == 1.0 after normalization
TEST(Normalize, FrobeniusNormAfterNormalization) {
    vnl_matrix<double> x(4, 2);
    x(0, 0)=1; x(0, 1)=2;
    x(1, 0)=3; x(1, 1)=4;
    x(2, 0)=5; x(2, 1)=6;
    x(3, 0)=7; x(3, 1)=8;
    vnl_vector<double> centroid;
    double scale;
    Normalize(x, centroid, scale);
    EXPECT_NEAR(x.frobenius_norm(), std::sqrt(4.0), 1e-10);
}

// The returned centroid equals the column-wise mean of the original x
TEST(Normalize, ReturnsCorrectCentroid) {
    vnl_matrix<double> x(3, 2);
    x(0, 0)=2; x(0, 1)=4;
    x(1, 0)=4; x(1, 1)=6;
    x(2, 0)=6; x(2, 1)=8;
    vnl_vector<double> centroid;
    double scale;
    Normalize(x, centroid, scale);
    EXPECT_NEAR(centroid[0], 4.0, 1e-10);
    EXPECT_NEAR(centroid[1], 6.0, 1e-10);
}

// The returned scale equals Frobenius-norm of centered x / sqrt(n)
TEST(Normalize, ReturnsCorrectScale) {
    vnl_matrix<double> x(3, 2);
    x(0, 0)=0; x(0, 1)=0;
    x(1, 0)=3; x(1, 1)=0;
    x(2, 0)=0; x(2, 1)=4;
    // centroid = (1, 4/3)
    // centered: (-1,-4/3), (2,-4/3), (-1,8/3)
    // ||centered||_F^2 = 1+16/9 + 4+16/9 + 1+64/9 = 6 + 96/9 = 6 + 32/3
    vnl_matrix<double> x_ref = x;
    vnl_vector<double> centroid_ref(2);
    centroid_ref[0] = 1.0; centroid_ref[1] = 4.0 / 3.0;
    for (int i = 0; i < 3; ++i) x_ref.set_row(i, x_ref.get_row(i) - centroid_ref);
    double expected_scale = x_ref.frobenius_norm() / std::sqrt(3.0);

    vnl_vector<double> centroid;
    double scale;
    Normalize(x, centroid, scale);
    EXPECT_NEAR(scale, expected_scale, 1e-10);
}

// Calling Normalize on an empty matrix must not crash
TEST(Normalize, EmptyMatrixNoOp) {
    vnl_matrix<double> x(0, 2);
    vnl_vector<double> centroid;
    double scale = 0.0;
    EXPECT_NO_THROW(Normalize(x, centroid, scale));
}

// ── Denormalize ──────────────────────────────────────────────────────────────

// Denormalize must exactly invert Normalize (round-trip)
TEST(Denormalize, InvertsNormalize) {
    vnl_matrix<double> x_orig(4, 2);
    x_orig(0, 0)=1; x_orig(0, 1)=2;
    x_orig(1, 0)=3; x_orig(1, 1)=4;
    x_orig(2, 0)=5; x_orig(2, 1)=6;
    x_orig(3, 0)=7; x_orig(3, 1)=8;
    vnl_matrix<double> x = x_orig;
    vnl_vector<double> centroid;
    double scale;
    Normalize(x, centroid, scale);
    Denormalize(x, centroid, scale);
    for (int i = 0; i < 4; ++i)
        for (int d = 0; d < 2; ++d)
            EXPECT_NEAR(x(i, d), x_orig(i, d), 1e-10)
                << "mismatch at (" << i << ", " << d << ")";
}

// Calling Denormalize on an empty matrix must not crash
TEST(Denormalize, EmptyMatrixNoOp) {
    vnl_matrix<double> x(0, 2);
    vnl_vector<double> centroid(2, 0.0);
    EXPECT_NO_THROW(Denormalize(x, centroid, 1.0));
}
