#include <cmath>

#include <gtest/gtest.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>

#include "gmmreg_rigid.h"
#include "gmmreg_rigid_func.h"

namespace {

// Small 2-D point sets used across multiple tests.
vnl_matrix<double> MakeModel2D() {
  vnl_matrix<double> m(3, 2);
  m(0, 0) = 0.0; m(0, 1) = 0.0;
  m(1, 0) = 1.0; m(1, 1) = 0.0;
  m(2, 0) = 0.5; m(2, 1) = 1.0;
  return m;
}

vnl_matrix<double> MakeScene2D() {
  vnl_matrix<double> s(4, 2);
  s(0, 0) = 0.1; s(0, 1) = 0.1;
  s(1, 0) = 0.9; s(1, 1) = 0.1;
  s(2, 0) = 0.5; s(2, 1) = 0.9;
  s(3, 0) = 0.5; s(3, 1) = 0.5;
  return s;
}

// Fixture: 2-D RigidRegistration loaded with synthetic point sets.
class RigidFunc2DTest : public ::testing::Test {
 protected:
  void SetUp() override {
    reg_.SetModelAndScene(MakeModel2D(), MakeScene2D());
    reg_.BuildTrees();
    func_.SetBase(&reg_);
    func_.SetScale(0.5);
  }
  gmmreg::RigidRegistration reg_;
  gmmreg::RigidFunc func_;
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
  EXPECT_DOUBLE_EQ(func.Eval(0.0, &g), 0.0);
  EXPECT_DOUBLE_EQ(g(0, 0), 0.0);
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
  gmmreg::RigidRegistration reg;
  reg.SetModelAndScene(MakeModel2D(), MakeScene2D());
  reg.BuildTrees();

  gmmreg::RigidFunc func;
  func.SetBase(&reg);
  EXPECT_EQ(func.get_number_of_unknowns(), 3);
}

TEST(RigidFuncSetBase, ThreeDimensionalGivesSevenUnknowns) {
  vnl_matrix<double> model(3, 3), scene(4, 3);
  model.fill(0.0); scene.fill(0.0);

  gmmreg::RigidRegistration reg;
  reg.SetModelAndScene(model, scene);
  reg.BuildTrees();

  gmmreg::RigidFunc func;
  func.SetBase(&reg);
  EXPECT_EQ(func.get_number_of_unknowns(), 7);
}

// ── SetModelAndScene ──────────────────────────────────────────────────────────

TEST(BaseSetModelAndScene, DimensionMismatchReturnsError) {
  gmmreg::RigidRegistration reg;
  vnl_matrix<double> model(3, 2), scene(4, 3);  // 2-D vs 3-D
  model.fill(0.0); scene.fill(0.0);
  EXPECT_LT(reg.SetModelAndScene(model, scene), 0);
}

TEST(BaseSetModelAndScene, SetsCorrectDimensions) {
  gmmreg::RigidRegistration reg;
  ASSERT_EQ(reg.SetModelAndScene(MakeModel2D(), MakeScene2D()), 0);
  // Verify via RigidFunc which reads m_ and d_ from Base.
  gmmreg::RigidFunc func;
  func.SetBase(&reg);
  EXPECT_EQ(func.get_number_of_unknowns(), 3);
}

// ── f ─────────────────────────────────────────────────────────────────────────

TEST_F(RigidFunc2DTest, IdentityReturnsFiniteNegative) {
  vnl_vector<double> x(3, 0.0);
  double val = func_.f(x);
  EXPECT_TRUE(std::isfinite(val));
  // Gauss cross-term > 0; Eval negates it.
  EXPECT_LT(val, 0.0);
}

TEST_F(RigidFunc2DTest, LargeTranslationNearZero) {
  vnl_vector<double> x_id(3, 0.0);
  double val_id = func_.f(x_id);

  // Move model 100 units away: cross-term ≈ 0 → f ≈ 0.
  vnl_vector<double> x_far(3, 0.0);
  x_far[0] = 100.0;
  double val_far = func_.f(x_far);

  EXPECT_LT(val_id, val_far);
}

// ── gradf ─────────────────────────────────────────────────────────────────────

TEST_F(RigidFunc2DTest, GradientAtIdentityIsFinite) {
  vnl_vector<double> x(3, 0.0);
  func_.f(x);
  vnl_vector<double> g(3, 0.0);
  func_.gradf(x, g);
  for (int i = 0; i < 3; ++i)
    EXPECT_TRUE(std::isfinite(g[i])) << "gradient[" << i << "] is not finite";
}

TEST_F(RigidFunc2DTest, MatchesNumericalGradient) {
  vnl_vector<double> x(3);
  x[0] = 0.1; x[1] = -0.1; x[2] = 0.05;

  // Analytical gradient (gradf uses gradient_ set by the preceding f call).
  func_.f(x);
  vnl_vector<double> ag(3, 0.0);
  func_.gradf(x, ag);

  // Numerical gradient via central differences.
  const double h = 1e-5;
  vnl_vector<double> ng(3, 0.0);
  for (int i = 0; i < 3; ++i) {
    vnl_vector<double> xp = x, xm = x;
    xp[i] += h; xm[i] -= h;
    ng[i] = (func_.f(xp) - func_.f(xm)) / (2.0 * h);
  }

  for (int i = 0; i < 3; ++i)
    EXPECT_NEAR(ag[i], ng[i], 1e-4) << "gradient mismatch at index " << i;
}
