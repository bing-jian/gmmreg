#include <memory>

#include <gtest/gtest.h>
#include <vnl/vnl_matrix.h>

#include "fgt_utils.h"
#include "gauss_transform.h"

using namespace gmmreg;

template <typename T>
std::unique_ptr<NanoflannTree<T>> MakeTree(const vnl_matrix<T>& pts) {
    auto tree = std::make_unique<NanoflannTree<T>>(pts);
    tree->tree.buildIndex();
    return tree;
}

// ── NanoflannTree ────────────────────────────────────────────────────────────

TEST(NanoflannTree, PointCount) {
    vnl_matrix<double> pts(5, 3);
    pts.fill(0.0);
    auto tree = MakeTree(pts);
    EXPECT_EQ(tree->matrix_adaptor.kdtree_get_point_count(), 5u);
}

// ── FastGaussTransform ───────────────────────────────────────────────────────

// Single identical point: all neighbors found, result ≈ 1.0
TEST(FastGaussTransform, IdenticalSinglePoint) {
    vnl_matrix<double> A(1, 2), B(1, 2);
    A(0, 0) = 1.0; A(0, 1) = 2.0;
    B = A;
    auto tree = MakeTree(B);
    vnl_matrix<double> grad(1, 2, 0.0);
    EXPECT_NEAR(FastGaussTransform(*tree, A, 1.0, grad), 1.0, 1e-4);
}

// Large scale captures all neighbors: FGT should closely match exact
TEST(FastGaussTransform, MatchesExactForLargeScale) {
    vnl_matrix<double> A(4, 2), B(4, 2);
    A(0,0)=0.1; A(0,1)=0.2;  A(1,0)=0.5; A(1,1)=0.3;
    A(2,0)=0.8; A(2,1)=0.7;  A(3,0)=0.3; A(3,1)=0.9;
    B(0,0)=0.2; B(0,1)=0.4;  B(1,0)=0.6; B(1,1)=0.1;
    B(2,0)=0.3; B(2,1)=0.8;  B(3,0)=0.7; B(3,1)=0.5;
    double scale = 5.0;

    auto tree = MakeTree(B);
    vnl_matrix<double> fast_grad(4, 2, 0.0), exact_grad(4, 2, 0.0);
    double fast  = FastGaussTransform(*tree, A, scale, fast_grad);
    double exact = GaussTransform(A, B, scale, exact_grad);

    EXPECT_NEAR(fast, exact, 1e-3);
    for (int i = 0; i < 4; ++i)
        for (int d = 0; d < 2; ++d)
            EXPECT_NEAR(fast_grad(i, d), exact_grad(i, d), 1e-3)
                << "Gradient mismatch at (" << i << ", " << d << ")";
}

// Result is always in (0, 1]
TEST(FastGaussTransform, ResultInRange) {
    vnl_matrix<double> A(3, 2), B(3, 2);
    A(0,0)=0.1; A(0,1)=0.2;  A(1,0)=0.5; A(1,1)=0.3;  A(2,0)=0.8; A(2,1)=0.7;
    B(0,0)=0.2; B(0,1)=0.4;  B(1,0)=0.6; B(1,1)=0.1;  B(2,0)=0.3; B(2,1)=0.8;
    auto tree = MakeTree(B);
    vnl_matrix<double> grad(3, 2, 0.0);
    double result = FastGaussTransform(*tree, A, 0.5, grad);
    EXPECT_GT(result, 0.0);
    EXPECT_LE(result, 1.0);
}

// Gradient matches central finite differences (large scale minimises cutoff effects)
TEST(FastGaussTransform, GradientNumericalCheck) {
    vnl_matrix<double> A(3, 2), B(4, 2);
    A(0,0)=0.1; A(0,1)=0.2;  A(1,0)=0.5; A(1,1)=0.3;  A(2,0)=0.8; A(2,1)=0.7;
    B(0,0)=0.2; B(0,1)=0.4;  B(1,0)=0.6; B(1,1)=0.1;
    B(2,0)=0.3; B(2,1)=0.9;  B(3,0)=0.7; B(3,1)=0.5;
    double scale = 5.0;

    auto tree = MakeTree(B);
    vnl_matrix<double> grad(3, 2, 0.0);
    FastGaussTransform(*tree, A, scale, grad);

    const double eps = 1e-5;
    for (int i = 0; i < 3; ++i) {
        for (int d = 0; d < 2; ++d) {
            vnl_matrix<double> A_plus = A, A_minus = A;
            A_plus(i, d) += eps;
            A_minus(i, d) -= eps;
            vnl_matrix<double> g(3, 2, 0.0);
            double fp = FastGaussTransform(*tree, A_plus, scale, g);
            g.fill(0.0);
            double fm = FastGaussTransform(*tree, A_minus, scale, g);
            double fd = (fp - fm) / (2.0 * eps);
            EXPECT_NEAR(grad(i, d), fd, 1e-4)
                << "Gradient mismatch at (" << i << ", " << d << ")";
        }
    }
}

// ── FastNeighborSearch ───────────────────────────────────────────────────────

// Large scale: every point finds all others, yielding m*m edges
TEST(FastNeighborSearch, LargeScaleAllPairsFound) {
    vnl_matrix<double> A(3, 2);
    A(0,0)=0.0; A(0,1)=0.0;
    A(1,0)=0.1; A(1,1)=0.1;
    A(2,0)=0.2; A(2,1)=0.2;
    auto tree = MakeTree(A);
    std::vector<std::pair<int, int>> edges;
    FastNeighborSearch(*tree, A, 100.0, &edges);
    EXPECT_EQ(edges.size(), 9u);
}

// Tiny scale with distant points: only self-pairs (dist=0) are within radius
TEST(FastNeighborSearch, SmallScaleOnlySelfEdges) {
    vnl_matrix<double> A(2, 2);
    A(0,0)=0.0;   A(0,1)=0.0;
    A(1,0)=100.0; A(1,1)=100.0;
    auto tree = MakeTree(A);
    std::vector<std::pair<int, int>> edges;
    FastNeighborSearch(*tree, A, 0.001, &edges);
    for (const auto& e : edges)
        EXPECT_EQ(e.first, e.second) << "Non-self edge found with tiny scale";
}

// ── FastSelfGaussTransform ───────────────────────────────────────────────────

// Large scale: self-transform via edges should match GaussTransform(A, A)
TEST(FastSelfGaussTransform, MatchesGaussTransformSelf) {
    vnl_matrix<double> A(4, 2);
    A(0,0)=0.1; A(0,1)=0.2;  A(1,0)=0.5; A(1,1)=0.3;
    A(2,0)=0.8; A(2,1)=0.7;  A(3,0)=0.3; A(3,1)=0.9;
    double scale = 5.0;

    auto tree = MakeTree(A);
    std::vector<std::pair<int, int>> edges;
    FastNeighborSearch(*tree, A, scale, &edges);

    vnl_matrix<double> fast_grad(4, 2, 0.0);
    double fast  = FastSelfGaussTransform(A, edges, scale, fast_grad);

    vnl_matrix<double> exact_grad(4, 2, 0.0);
    double exact = GaussTransform(A, A, scale, exact_grad);

    EXPECT_NEAR(fast, exact, 1e-3);
}

// Gradient matches central finite differences (edge list held fixed during perturbation)
TEST(FastSelfGaussTransform, GradientNumericalCheck) {
    vnl_matrix<double> A(4, 2);
    A(0,0)=0.1; A(0,1)=0.2;  A(1,0)=0.5; A(1,1)=0.3;
    A(2,0)=0.8; A(2,1)=0.7;  A(3,0)=0.3; A(3,1)=0.9;
    double scale = 5.0;

    auto tree = MakeTree(A);
    std::vector<std::pair<int, int>> edges;
    FastNeighborSearch(*tree, A, scale, &edges);

    vnl_matrix<double> grad(4, 2, 0.0);
    FastSelfGaussTransform(A, edges, scale, grad);

    const double eps = 1e-5;
    for (int i = 0; i < 4; ++i) {
        for (int d = 0; d < 2; ++d) {
            vnl_matrix<double> A_plus = A, A_minus = A;
            A_plus(i, d) += eps;
            A_minus(i, d) -= eps;
            vnl_matrix<double> g(4, 2, 0.0);
            double fp = FastSelfGaussTransform(A_plus,  edges, scale, g);
            g.fill(0.0);
            double fm = FastSelfGaussTransform(A_minus, edges, scale, g);
            double fd = (fp - fm) / (2.0 * eps);
            EXPECT_NEAR(grad(i, d), fd, 1e-4)
                << "Gradient mismatch at (" << i << ", " << d << ")";
        }
    }
}
