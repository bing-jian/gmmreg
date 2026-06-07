#include <vector>

#include <gtest/gtest.h>
#include <vnl/vnl_matrix.h>

#include "match_utils.h"

using namespace gmmreg;

// ── ComputeSquaredDistanceMatrix ─────────────────────────────────────────────

// Output size is m×n
TEST(ComputeSquaredDistanceMatrix, OutputSize) {
    vnl_matrix<double> A(3, 2), B(4, 2);
    A.fill(0); B.fill(1);
    vnl_matrix<double> D;
    ComputeSquaredDistanceMatrix(A, B, D);
    EXPECT_EQ(D.rows(), 3u);
    EXPECT_EQ(D.cols(), 4u);
}

// When A == B the diagonal entries are 0
TEST(ComputeSquaredDistanceMatrix, SelfDistanceIsZero) {
    vnl_matrix<double> A(3, 2);
    A(0,0)=1; A(0,1)=2;
    A(1,0)=3; A(1,1)=4;
    A(2,0)=5; A(2,1)=6;
    vnl_matrix<double> D;
    ComputeSquaredDistanceMatrix(A, A, D);
    for (int i = 0; i < 3; ++i)
        EXPECT_NEAR(D(i, i), 0.0, 1e-10);
}

// Known values: A=[0,0], B=[3,4] → D[0,0] = 9+16 = 25
TEST(ComputeSquaredDistanceMatrix, KnownValue) {
    vnl_matrix<double> A(1, 2), B(1, 2);
    A(0,0)=0; A(0,1)=0;
    B(0,0)=3; B(0,1)=4;
    vnl_matrix<double> D;
    ComputeSquaredDistanceMatrix(A, B, D);
    EXPECT_NEAR(D(0, 0), 25.0, 1e-10);
}

// D_AB[i,j] == D_BA[j,i] (squared distance is symmetric)
TEST(ComputeSquaredDistanceMatrix, Symmetry) {
    vnl_matrix<double> A(2, 2), B(3, 2);
    A(0,0)=0; A(0,1)=0; A(1,0)=1; A(1,1)=1;
    B(0,0)=2; B(0,1)=0; B(1,0)=0; B(1,1)=3; B(2,0)=1; B(2,1)=2;
    vnl_matrix<double> DAB, DBA;
    ComputeSquaredDistanceMatrix(A, B, DAB);
    ComputeSquaredDistanceMatrix(B, A, DBA);
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 3; ++j)
            EXPECT_NEAR(DAB(i, j), DBA(j, i), 1e-10);
}

// ── SelectPoints ─────────────────────────────────────────────────────────────

// Returns the correct rows and the right count
TEST(SelectPoints, PicksCorrectRows) {
    vnl_matrix<double> pts(4, 2);
    pts(0,0)=10; pts(0,1)=11;
    pts(1,0)=20; pts(1,1)=21;
    pts(2,0)=30; pts(2,1)=31;
    pts(3,0)=40; pts(3,1)=41;
    std::vector<int> idx = {0, 2, 3};
    vnl_matrix<double> sel;
    int n = SelectPoints(pts, idx, sel);
    EXPECT_EQ(n, 3);
    EXPECT_EQ(sel.rows(), 3u);
    EXPECT_NEAR(sel(0,0), 10.0, 1e-10);
    EXPECT_NEAR(sel(1,0), 30.0, 1e-10);
    EXPECT_NEAR(sel(2,0), 40.0, 1e-10);
}

// Empty index vector → empty selected, return 0
TEST(SelectPoints, EmptyIndex) {
    vnl_matrix<double> pts(3, 2);
    pts.fill(1.0);
    std::vector<int> idx;
    vnl_matrix<double> sel;
    int n = SelectPoints(pts, idx, sel);
    EXPECT_EQ(n, 0);
    EXPECT_EQ(sel.rows(), 0u);
}

// ── FindNearestPairs ─────────────────────────────────────────────────────────

// Diagonal-dominant matrix → diagonal pairs selected
TEST(FindNearestPairs, DiagonalMatch) {
    vnl_matrix<double> dist(2, 2);
    dist(0,0)=0.1; dist(0,1)=5.0;
    dist(1,0)=5.0; dist(1,1)=0.1;
    vnl_matrix<int> pairs;
    FindNearestPairs(dist, pairs, 1.0);
    EXPECT_EQ(pairs.cols(), 2u);
    // pairs row 0 = row indices, row 1 = col indices
    bool got_00 = (pairs(0,0)==0 && pairs(1,0)==0);
    bool got_11 = (pairs(0,1)==1 && pairs(1,1)==1);
    EXPECT_TRUE(got_00);
    EXPECT_TRUE(got_11);
}

// Threshold below all distances → no pairs found
TEST(FindNearestPairs, ThresholdTooLow) {
    vnl_matrix<double> dist(2, 2);
    dist(0,0)=1.0; dist(0,1)=2.0;
    dist(1,0)=3.0; dist(1,1)=4.0;
    vnl_matrix<int> pairs;
    FindNearestPairs(dist, pairs, 0.5);
    EXPECT_EQ(pairs.cols(), 0u);
}

// Each column appears at most once in the output
TEST(FindNearestPairs, ColumnSelectedAtMostOnce) {
    // Both rows have their minimum in column 0
    vnl_matrix<double> dist(2, 2);
    dist(0,0)=0.1; dist(0,1)=5.0;
    dist(1,0)=0.2; dist(1,1)=5.0;
    vnl_matrix<int> pairs;
    FindNearestPairs(dist, pairs, 1.0);
    // Row 0 grabs col 0; row 1 cannot grab col 0 again → only 1 pair
    EXPECT_EQ(pairs.cols(), 1u);
}

// ── PickIndices ──────────────────────────────────────────────────────────────

// Diagonal-dominant matrix → diagonal indices selected
TEST(PickIndices, DiagonalMatch) {
    vnl_matrix<double> dist(2, 2);
    dist(0,0)=0.1; dist(0,1)=5.0;
    dist(1,0)=5.0; dist(1,1)=0.1;
    std::vector<int> row_idx, col_idx;
    PickIndices(dist, row_idx, col_idx, 1.0);
    EXPECT_EQ(row_idx.size(), 2u);
    EXPECT_EQ(col_idx.size(), 2u);
    EXPECT_EQ(row_idx[0], 0); EXPECT_EQ(col_idx[0], 0);
    EXPECT_EQ(row_idx[1], 1); EXPECT_EQ(col_idx[1], 1);
}

// Threshold below all values → no indices
TEST(PickIndices, ThresholdTooLow) {
    vnl_matrix<double> dist(2, 2);
    dist.fill(1.0);
    std::vector<int> row_idx, col_idx;
    PickIndices(dist, row_idx, col_idx, 0.5);
    EXPECT_TRUE(row_idx.empty());
    EXPECT_TRUE(col_idx.empty());
}

// ── FindWorkingPairs ─────────────────────────────────────────────────────────

// Transformed_M == S → all points pair up (threshold > 0)
TEST(FindWorkingPairs, ExactMatchReturnsAllPairs) {
    vnl_matrix<double> M(3, 2), S(3, 2), TM(3, 2);
    M(0,0)=0; M(0,1)=0; M(1,0)=1; M(1,1)=0; M(2,0)=2; M(2,1)=0;
    S(0,0)=0; S(0,1)=1; S(1,0)=1; S(1,1)=1; S(2,0)=2; S(2,1)=1;
    TM = S;  // perfectly aligned
    vnl_matrix<double> wM, wS;
    int n = FindWorkingPairs(M, S, TM, 0.01, wM, wS);
    EXPECT_EQ(n, 3);
    EXPECT_EQ(wM.rows(), 3u);
    EXPECT_EQ(wS.rows(), 3u);
}

// Transformed_M far from S → no pairs below threshold
TEST(FindWorkingPairs, NoMatchWhenFarApart) {
    vnl_matrix<double> M(2, 2), S(2, 2), TM(2, 2);
    M.fill(0); S.fill(0); TM.fill(100.0);
    vnl_matrix<double> wM, wS;
    int n = FindWorkingPairs(M, S, TM, 1.0, wM, wS);
    EXPECT_EQ(n, 0);
}
