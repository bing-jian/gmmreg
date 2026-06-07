#include <cmath>

#include <gtest/gtest.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>

#include "rbf_utils.h"

using namespace gmmreg;

// ── GaussianAffinityMatrix ───────────────────────────────────────────────────

// Identical points → affinity = 1.0
TEST(GaussianAffinityMatrix, IdenticalPoints) {
    double A[] = {1.0, 2.0};
    double B[] = {1.0, 2.0};
    double dist[1];
    GaussianAffinityMatrix(A, B, 1, 1, 2, 1.0, dist);
    EXPECT_NEAR(dist[0], 1.0, 1e-10);
}

// Very distant points with small scale → affinity ≈ 0
TEST(GaussianAffinityMatrix, DistantPointsSmallScale) {
    double A[] = {0.0, 0.0};
    double B[] = {1000.0, 1000.0};
    double dist[1];
    GaussianAffinityMatrix(A, B, 1, 1, 2, 0.01, dist);
    EXPECT_NEAR(dist[0], 0.0, 1e-10);
}

// Larger scale → larger affinity for the same pair of points
TEST(GaussianAffinityMatrix, LargerScaleHigherAffinity) {
    double A[] = {0.0, 0.0};
    double B[] = {1.0, 0.0};
    double d_tight[1], d_wide[1];
    GaussianAffinityMatrix(A, B, 1, 1, 2, 0.1, d_tight);
    GaussianAffinityMatrix(A, B, 1, 1, 2, 10.0, d_wide);
    EXPECT_LT(d_tight[0], d_wide[0]);
}

// Swapping A and B transposes the output matrix
TEST(GaussianAffinityMatrix, SwapTransposesResult) {
    double A[] = {0.0, 0.0, 1.0, 0.0};
    double B[] = {0.5, 0.5, 1.5, 0.5};
    double dAB[4], dBA[4];
    GaussianAffinityMatrix(A, B, 2, 2, 2, 1.0, dAB);
    GaussianAffinityMatrix(B, A, 2, 2, 2, 1.0, dBA);
    // dAB[i*2+j] == dBA[j*2+i]
    EXPECT_NEAR(dAB[0], dBA[0], 1e-10);
    EXPECT_NEAR(dAB[1], dBA[2], 1e-10);
    EXPECT_NEAR(dAB[2], dBA[1], 1e-10);
    EXPECT_NEAR(dAB[3], dBA[3], 1e-10);
}

// ── ComputeTPSKernel ─────────────────────────────────────────────────────────

// U is m×n and K is n×n
TEST(ComputeTPSKernel, OutputSizes) {
    vnl_matrix<double> model(4, 2), ctrl(3, 2);
    model.fill(0); ctrl.fill(0);
    // Scatter points so they are distinct and r > eps
    for (int i = 0; i < 4; ++i) { model(i, 0) = i * 1.0; model(i, 1) = 0.0; }
    for (int i = 0; i < 3; ++i) { ctrl(i, 0) = i * 1.5; ctrl(i, 1) = 1.0; }
    vnl_matrix<double> U, K;
    ComputeTPSKernel(model, ctrl, U, K);
    EXPECT_EQ(U.rows(), 4u);
    EXPECT_EQ(U.cols(), 3u);
    EXPECT_EQ(K.rows(), 3u);
    EXPECT_EQ(K.cols(), 3u);
}

// Diagonal of K is 0 (same point: r = 0, skipped by eps guard)
TEST(ComputeTPSKernel, KDiagonalZero2D) {
    vnl_matrix<double> pts(3, 2);
    pts(0, 0) = 0.0; pts(0, 1) = 0.0;
    pts(1, 0) = 1.0; pts(1, 1) = 0.0;
    pts(2, 0) = 0.0; pts(2, 1) = 1.0;
    vnl_matrix<double> U, K;
    ComputeTPSKernel(pts, pts, U, K);
    for (int i = 0; i < 3; ++i)
        EXPECT_NEAR(K(i, i), 0.0, 1e-10) << "K(" << i << "," << i << ") != 0";
}

// K must be symmetric
TEST(ComputeTPSKernel, KSymmetric2D) {
    vnl_matrix<double> pts(3, 2);
    pts(0, 0) = 0.0; pts(0, 1) = 0.0;
    pts(1, 0) = 1.0; pts(1, 1) = 0.0;
    pts(2, 0) = 0.0; pts(2, 1) = 1.0;
    vnl_matrix<double> U, K;
    ComputeTPSKernel(pts, pts, U, K);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            EXPECT_NEAR(K(i, j), K(j, i), 1e-10);
}

// ── ComputeGaussianKernel ────────────────────────────────────────────────────

// When model == ctrl_pts the function sets K = G
TEST(ComputeGaussianKernel, GEqualsKWhenModelIsCtrlPts) {
    vnl_matrix<double> pts(3, 2);
    pts(0, 0) = 0.0; pts(0, 1) = 0.0;
    pts(1, 0) = 1.0; pts(1, 1) = 0.0;
    pts(2, 0) = 0.5; pts(2, 1) = 0.5;
    vnl_matrix<double> G, K;
    ComputeGaussianKernel(pts, pts, G, K, 1.0);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            EXPECT_NEAR(G(i, j), K(i, j), 1e-10);
}

// Diagonal of G is 1.0 when model == ctrl_pts (distance 0 → exp(0) = 1)
TEST(ComputeGaussianKernel, DiagonalIsOneWhenModelIsCtrlPts) {
    vnl_matrix<double> pts(3, 2);
    pts(0, 0) = 0.0; pts(0, 1) = 0.0;
    pts(1, 0) = 1.0; pts(1, 1) = 0.0;
    pts(2, 0) = 0.5; pts(2, 1) = 0.5;
    vnl_matrix<double> G, K;
    ComputeGaussianKernel(pts, pts, G, K, 1.0);
    for (int i = 0; i < 3; ++i)
        EXPECT_NEAR(G(i, i), 1.0, 1e-10) << "G(" << i << "," << i << ") != 1";
}

// K is symmetric (ctrl_pts self-affinity)
TEST(ComputeGaussianKernel, KSymmetric) {
    vnl_matrix<double> model(4, 2), ctrl(3, 2);
    model(0, 0)=0.1; model(0, 1)=0.2;
    model(1, 0)=0.5; model(1, 1)=0.3;
    model(2, 0)=0.8; model(2, 1)=0.7;
    model(3, 0)=0.3; model(3, 1)=0.9;
    ctrl(0, 0)=0.2; ctrl(0, 1)=0.4;
    ctrl(1, 0)=0.6; ctrl(1, 1)=0.1;
    ctrl(2, 0)=0.3; ctrl(2, 1)=0.8;
    vnl_matrix<double> G, K;
    ComputeGaussianKernel(model, ctrl, G, K, 1.0);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            EXPECT_NEAR(K(i, j), K(j, i), 1e-10);
}
