#include <cmath>

#include <gtest/gtest.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>

#include "gmmreg_base.h"
#include "gmmreg_rigid_func.h"

namespace {

// Minimal 2-D Base: implements the same rigid transform as RigidRegistration.
class TestBase2D : public gmmreg::Base {
 public:
  TestBase2D() {
    d_ = 2; m_ = 3; s_ = 4;

    model_.set_size(3, 2);
    model_(0, 0) = 0.0; model_(0, 1) = 0.0;
    model_(1, 0) = 1.0; model_(1, 1) = 0.0;
    model_(2, 0) = 0.5; model_(2, 1) = 1.0;

    scene_.set_size(4, 2);
    scene_(0, 0) = 0.1; scene_(0, 1) = 0.1;
    scene_(1, 0) = 0.9; scene_(1, 1) = 0.1;
    scene_(2, 0) = 0.5; scene_(2, 1) = 0.9;
    scene_(3, 0) = 0.5; scene_(3, 1) = 0.5;

    transformed_model_ = model_;

#ifdef USE_KDTREE
    scene_tree_ = std::make_unique<gmmreg::NanoflannTree<double>>(scene_);
    scene_tree_->tree.buildIndex();
#endif
  }

  // Mirrors ConvertParamToRotationAndTranslation (2-D):
  //   transformed_model_ = model_ * R.T + ones * [tx, ty]
  void PerformTransform(const vnl_vector<double>& x) override {
    double tx = x[0], ty = x[1], theta = x[2];
    double c = std::cos(theta), s = std::sin(theta);
    // R.T = [[c, s], [-s, c]]
    vnl_matrix<double> RT(2, 2);
    RT(0, 0) = c;  RT(0, 1) = s;
    RT(1, 0) = -s; RT(1, 1) = c;
    transformed_model_ = model_ * RT;
    for (int i = 0; i < m_; ++i) {
      transformed_model_(i, 0) += tx;
      transformed_model_(i, 1) += ty;
    }
  }

  double BendingEnergy() override { return 0.0; }
  void ComputeGradient(const double, const vnl_matrix<double>&,
                       vnl_matrix<double>&) override {}

 private:
  void PrepareOwnOptions(const char*) override {}
  void PrepareBasisKernel() override {}
  int SetInitParams(const char*) override { return 0; }
  void SaveResults(const char*, const vnl_vector<double>&) override {}
  void StartRegistration(vnl_vector<double>&) override {}
};

// Minimal 3-D Base: only used to verify SetBase gives 7 unknowns.
class TestBase3D : public gmmreg::Base {
 public:
  TestBase3D() {
    d_ = 3; m_ = 3; s_ = 3;

    model_.set_size(3, 3);
    model_(0, 0) = 1.0; model_(0, 1) = 0.0; model_(0, 2) = 0.0;
    model_(1, 0) = 0.0; model_(1, 1) = 1.0; model_(1, 2) = 0.0;
    model_(2, 0) = 0.0; model_(2, 1) = 0.0; model_(2, 2) = 1.0;

    scene_ = model_;
    transformed_model_ = model_;

#ifdef USE_KDTREE
    scene_tree_ = std::make_unique<gmmreg::NanoflannTree<double>>(scene_);
    scene_tree_->tree.buildIndex();
#endif
  }

  void PerformTransform(const vnl_vector<double>&) override {}
  double BendingEnergy() override { return 0.0; }
  void ComputeGradient(const double, const vnl_matrix<double>&,
                       vnl_matrix<double>&) override {}

 private:
  void PrepareOwnOptions(const char*) override {}
  void PrepareBasisKernel() override {}
  int SetInitParams(const char*) override { return 0; }
  void SaveResults(const char*, const vnl_vector<double>&) override {}
  void StartRegistration(vnl_vector<double>&) override {}
};

}  // namespace

// ── Eval ─────────────────────────────────────────────────────────────────────

TEST(RigidFuncEval, NegatesEnergyAndGradient) {
  gmmreg::RigidFunc func;
  vnl_matrix<double> g(3, 2);
  g(0, 0) = 1.0; g(0, 1) = 2.0;
  g(1, 0) = 3.0; g(1, 1) = 4.0;
  g(2, 0) = 5.0; g(2, 1) = 6.0;

  double result = func.Eval(7.5, &g);
  EXPECT_DOUBLE_EQ(result, -7.5);
  EXPECT_DOUBLE_EQ(g(0, 0), -1.0);
  EXPECT_DOUBLE_EQ(g(1, 1), -4.0);
  EXPECT_DOUBLE_EQ(g(2, 1), -6.0);
}

TEST(RigidFuncEval, ZeroEnergyZeroGradient) {
  gmmreg::RigidFunc func;
  vnl_matrix<double> g(2, 2);
  g.fill(0.0);
  double result = func.Eval(0.0, &g);
  EXPECT_DOUBLE_EQ(result, 0.0);
  EXPECT_DOUBLE_EQ(g(0, 0), 0.0);
  EXPECT_DOUBLE_EQ(g(1, 1), 0.0);
}

// ── Scale ─────────────────────────────────────────────────────────────────────

TEST(RigidFuncScale, GetAfterSet) {
  gmmreg::RigidFunc func;
  func.SetScale(1.5);
  EXPECT_DOUBLE_EQ(func.GetScale(), 1.5);
  func.SetScale(0.01);
  EXPECT_DOUBLE_EQ(func.GetScale(), 0.01);
}

// ── SetBase ───────────────────────────────────────────────────────────────────

TEST(RigidFuncSetBase, TwoDimensionalGivesThreeUnknowns) {
  TestBase2D base;
  gmmreg::RigidFunc func;
  func.SetBase(&base);
  EXPECT_EQ(func.get_number_of_unknowns(), 3);
}

TEST(RigidFuncSetBase, ThreeDimensionalGivesSevenUnknowns) {
  TestBase3D base;
  gmmreg::RigidFunc func;
  func.SetBase(&base);
  EXPECT_EQ(func.get_number_of_unknowns(), 7);
}

// ── f ─────────────────────────────────────────────────────────────────────────

TEST(RigidFuncF, IdentityReturnsFiniteNegative) {
  TestBase2D base;
  gmmreg::RigidFunc func;
  func.SetBase(&base);
  func.SetScale(0.5);

  vnl_vector<double> x(3, 0.0);
  double val = func.f(x);
  EXPECT_TRUE(std::isfinite(val));
  // Gauss transform cross-term > 0; Eval negates it → f < 0
  EXPECT_LT(val, 0.0);
}

TEST(RigidFuncF, LargeTranslationNearZero) {
  TestBase2D base;
  gmmreg::RigidFunc func;
  func.SetBase(&base);
  func.SetScale(0.5);

  // Identity: model overlaps scene → large negative value
  vnl_vector<double> x_id(3, 0.0);
  double val_id = func.f(x_id);

  // Move model 100 units away: cross-term ≈ 0 → f ≈ 0
  vnl_vector<double> x_far(3, 0.0);
  x_far[0] = 100.0;
  double val_far = func.f(x_far);

  EXPECT_LT(val_id, val_far);
}

// ── gradf ─────────────────────────────────────────────────────────────────────

TEST(RigidFuncGradf, MatchesNumericalGradient2D) {
  TestBase2D base;
  gmmreg::RigidFunc func;
  func.SetBase(&base);
  func.SetScale(0.5);

  vnl_vector<double> x(3);
  x[0] = 0.1; x[1] = -0.1; x[2] = 0.05;

  // Analytical gradient: gradf uses gradient_ set by the preceding f() call.
  func.f(x);
  vnl_vector<double> ag(3, 0.0);
  func.gradf(x, ag);

  // Numerical gradient via central differences.
  const double h = 1e-5;
  vnl_vector<double> ng(3, 0.0);
  for (int i = 0; i < 3; ++i) {
    vnl_vector<double> xp = x, xm = x;
    xp[i] += h; xm[i] -= h;
    ng[i] = (func.f(xp) - func.f(xm)) / (2.0 * h);
  }

  for (int i = 0; i < 3; ++i)
    EXPECT_NEAR(ag[i], ng[i], 1e-4) << "gradient mismatch at index " << i;
}

TEST(RigidFuncGradf, GradientAtIdentityIsFinite) {
  TestBase2D base;
  gmmreg::RigidFunc func;
  func.SetBase(&base);
  func.SetScale(0.5);

  vnl_vector<double> x(3, 0.0);
  func.f(x);
  vnl_vector<double> g(3, 0.0);
  func.gradf(x, g);

  for (int i = 0; i < 3; ++i)
    EXPECT_TRUE(std::isfinite(g[i])) << "gradient[" << i << "] is not finite";
}
