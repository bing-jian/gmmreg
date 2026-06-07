#include <cmath>

#include <gtest/gtest.h>
#include <vnl/vnl_matrix.h>

#include "em_utils.h"

using namespace gmmreg;

// ── Column normalisation ─────────────────────────────────────────────────────

// Without outliers each column of P must sum exactly to 1.
TEST(ComputeP, ColumnSumsToOneNoOutliers) {
    vnl_matrix<double> x(3, 2), y(2, 2);
    x(0,0)=0.0; x(0,1)=0.0;
    x(1,0)=1.0; x(1,1)=0.0;
    x(2,0)=0.5; x(2,1)=0.5;
    y(0,0)=0.2; y(0,1)=0.3;
    y(1,0)=0.8; y(1,1)=0.7;
    vnl_matrix<double> P(3, 2, 0.0);
    double E = 0.0;

    ComputeP(x, y, P, E, 1.0, 0);

    for (int j = 0; j < 2; ++j) {
        double col_sum = 0.0;
        for (int i = 0; i < 3; ++i) col_sum += P(i, j);
        EXPECT_NEAR(col_sum, 1.0, 1e-10) << "Column " << j << " does not sum to 1";
    }
}

// With outliers the outlier term is added to each column denominator, so
// each column sum must be strictly less than 1.
TEST(ComputeP, ColumnSumsLessThanOneWithOutliers) {
    vnl_matrix<double> x(3, 2), y(2, 2);
    x(0,0)=0.0; x(0,1)=0.0;
    x(1,0)=1.0; x(1,1)=0.0;
    x(2,0)=0.5; x(2,1)=0.5;
    y(0,0)=0.2; y(0,1)=0.3;
    y(1,0)=0.8; y(1,1)=0.7;
    vnl_matrix<double> P(3, 2, 0.0);
    double E = 0.0;

    ComputeP(x, y, P, E, 1.0, 1);

    for (int j = 0; j < 2; ++j) {
        double col_sum = 0.0;
        for (int i = 0; i < 3; ++i) col_sum += P(i, j);
        EXPECT_LT(col_sum, 1.0) << "Column " << j << " should be < 1 with outliers";
    }
}

// ── Identical points ─────────────────────────────────────────────────────────

// Single identical point, no outliers: P must be 1.0.
TEST(ComputeP, SingleIdenticalPoint) {
    vnl_matrix<double> x(1, 2), y(1, 2);
    x(0,0)=1.0; x(0,1)=2.0;
    y = x;
    vnl_matrix<double> P(1, 1, 0.0);
    double E = 0.0;

    ComputeP(x, y, P, E, 1.0, 0);

    EXPECT_NEAR(P(0, 0), 1.0, 1e-10);
}

// When x == y (same point set) closer matches should dominate each column.
TEST(ComputeP, IdenticalSetsNearDiagonalDominance) {
    vnl_matrix<double> x(2, 2);
    x(0,0)=0.0; x(0,1)=0.0;
    x(1,0)=5.0; x(1,1)=0.0;
    vnl_matrix<double> y = x;
    vnl_matrix<double> P(2, 2, 0.0);
    double E = 0.0;

    ComputeP(x, y, P, E, 0.5, 0);

    // Each point should predominantly correspond to itself.
    EXPECT_GT(P(0, 0), P(1, 0));  // x[0] closer to y[0] than y[1]
    EXPECT_GT(P(1, 1), P(0, 1));  // x[1] closer to y[1] than y[0]
}

// ── Log-likelihood ───────────────────────────────────────────────────────────

// For a single pair of identical points E should be 0 (log(1) = 0).
TEST(ComputeP, LogLikelihoodIdenticalPoints) {
    vnl_matrix<double> x(1, 2), y(1, 2);
    x(0,0)=0.0; x(0,1)=0.0;
    y = x;
    vnl_matrix<double> P(1, 1, 0.0);
    double E = 0.0;

    ComputeP(x, y, P, E, 1.0, 0);

    EXPECT_NEAR(E, 0.0, 1e-10);
}

// E should equal -sum_j log(column_sum_j) for a known two-point case.
TEST(ComputeP, LogLikelihoodKnownValue) {
    // x has one point, y has one point at distance sqrt(2); sigma=1
    // P(0,0) = exp(-2 / (2*1)) = exp(-1)
    // column_sum = exp(-1), E = -log(exp(-1)) = 1.0
    vnl_matrix<double> x(1, 2), y(1, 2);
    x(0,0)=0.0; x(0,1)=0.0;
    y(0,0)=1.0; y(0,1)=1.0;
    vnl_matrix<double> P(1, 1, 0.0);
    double E = 0.0;

    ComputeP(x, y, P, E, 1.0, 0);

    EXPECT_NEAR(E, 1.0, 1e-10);
}

// ── Degenerate sigma (regression test for P.fill(0) fix) ────────────────────

// Extremely small sigma with distant points: column_sum falls below 1e-12,
// P must be zeroed out rather than left with unnormalised values.
TEST(ComputeP, DegenerateSigmaZerosP) {
    vnl_matrix<double> x(1, 2), y(1, 2);
    x(0,0)=0.0;   x(0,1)=0.0;
    y(0,0)=100.0; y(0,1)=100.0;
    vnl_matrix<double> P(1, 1, 0.0);
    double E = 0.0;

    ComputeP(x, y, P, E, 1e-10, 0);

    EXPECT_NEAR(P(0, 0), 0.0, 1e-30);
}

// ── Sigma effect ─────────────────────────────────────────────────────────────

// Larger sigma → more uniform correspondence (smaller max column entry).
TEST(ComputeP, LargerSigmaMoreUniform) {
    vnl_matrix<double> x(3, 2), y(3, 2);
    x(0,0)=0.0; x(0,1)=0.0;
    x(1,0)=1.0; x(1,1)=0.0;
    x(2,0)=2.0; x(2,1)=0.0;
    y = x;
    vnl_matrix<double> P_tight(3, 3, 0.0), P_wide(3, 3, 0.0);
    double E = 0.0;

    ComputeP(x, y, P_tight, E, 0.1, 0);
    ComputeP(x, y, P_wide,  E, 10.0, 0);

    double max_tight = P_tight.max_value();
    double max_wide  = P_wide.max_value();
    EXPECT_GT(max_tight, max_wide);
}
