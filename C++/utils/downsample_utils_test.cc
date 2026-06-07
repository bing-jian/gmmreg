#include <cmath>
#include <memory>

#include <gtest/gtest.h>
#include <vnl/vnl_matrix.h>

#include "downsample_utils.h"

using namespace gmmreg;

template <typename T>
std::unique_ptr<NanoflannTree<T>> MakeTree(const vnl_matrix<T>& pts) {
    auto tree = std::make_unique<NanoflannTree<T>>(pts);
    tree->tree.buildIndex();
    return tree;
}

// Output row count must not exceed input row count
TEST(CreateSparseNodes, OutputNotLargerThanInput) {
    vnl_matrix<double> pts(10, 2);
    for (int i = 0; i < 10; ++i) { pts(i, 0) = i * 1.0; pts(i, 1) = 0.0; }
    auto result = CreateSparseNodes(*MakeTree(pts), 0.5);
    EXPECT_LE(result.rows(), 10u);
    EXPECT_EQ(result.cols(), 2u);
}

// Radius larger than the whole point set → exactly one output point
TEST(CreateSparseNodes, LargeRadiusGivesOnePoint) {
    vnl_matrix<double> pts(5, 2);
    for (int i = 0; i < 5; ++i) { pts(i, 0) = i * 1.0; pts(i, 1) = 0.0; }
    // max spread is 4; radius 10 covers all points from any seed
    auto result = CreateSparseNodes(*MakeTree(pts), 10.0);
    EXPECT_EQ(result.rows(), 1u);
}

// Radius smaller than minimum inter-point spacing → all points kept
TEST(CreateSparseNodes, SmallRadiusKeepsAllPoints) {
    vnl_matrix<double> pts(5, 2);
    for (int i = 0; i < 5; ++i) { pts(i, 0) = i * 10.0; pts(i, 1) = 0.0; }
    // spacing = 10; radius 0.1 cannot remove any neighbour
    auto result = CreateSparseNodes(*MakeTree(pts), 0.1);
    EXPECT_EQ(result.rows(), 5u);
}

// Every row of the output must appear verbatim in the input
TEST(CreateSparseNodes, OutputIsSubsetOfInput) {
    vnl_matrix<double> pts(6, 2);
    for (int i = 0; i < 6; ++i) { pts(i, 0) = i * 2.0; pts(i, 1) = 0.0; }
    auto result = CreateSparseNodes(*MakeTree(pts), 1.0);
    for (int r = 0; r < result.rows(); ++r) {
        bool found = false;
        for (int i = 0; i < 6 && !found; ++i) {
            found = (std::abs(result(r, 0) - pts(i, 0)) < 1e-10 &&
                     std::abs(result(r, 1) - pts(i, 1)) < 1e-10);
        }
        EXPECT_TRUE(found) << "output row " << r << " not found in input";
    }
}

// No two output points may be closer than dist_radius to each other
TEST(CreateSparseNodes, OutputPointsMutuallyBeyondRadius) {
    vnl_matrix<double> pts(15, 2);
    for (int i = 0; i < 15; ++i) { pts(i, 0) = i * 1.0; pts(i, 1) = 0.0; }
    const double radius = 2.5;
    auto result = CreateSparseNodes(*MakeTree(pts), radius);
    int n = result.rows();
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            double dx = result(i, 0) - result(j, 0);
            double dy = result(i, 1) - result(j, 1);
            double dist = std::sqrt(dx*dx + dy*dy);
            EXPECT_GT(dist, radius - 1e-6)
                << "output points " << i << " and " << j
                << " are closer than radius: " << dist;
        }
    }
}
