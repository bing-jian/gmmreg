#include <cmath>

#include <gtest/gtest.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>

#include "transform_utils.h"

using namespace gmmreg;

static const double kEps = 1e-10;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static vnl_matrix<double> Identity44() {
  vnl_matrix<double> T(4, 4, 0.0);
  T(0, 0) = T(1, 1) = T(2, 2) = T(3, 3) = 1.0;
  return T;
}

static vnl_matrix<double> RotZ(double deg) {
  vnl_matrix<double> T = Identity44();
  double r = deg * M_PI / 180.0;
  T(0, 0) =  std::cos(r);  T(0, 1) = -std::sin(r);
  T(1, 0) =  std::sin(r);  T(1, 1) =  std::cos(r);
  return T;
}

static vnl_matrix<double> RotX(double deg) {
  vnl_matrix<double> T = Identity44();
  double r = deg * M_PI / 180.0;
  T(1, 1) =  std::cos(r);  T(1, 2) = -std::sin(r);
  T(2, 1) =  std::sin(r);  T(2, 2) =  std::cos(r);
  return T;
}

static vnl_matrix<double> Trans(double tx, double ty, double tz) {
  vnl_matrix<double> T = Identity44();
  T(0, 3) = tx;  T(1, 3) = ty;  T(2, 3) = tz;
  return T;
}

// ── RelativeTransform ─────────────────────────────────────────────────────────

TEST(RelativeTransform, BothIdentity) {
  vnl_matrix<double> I = Identity44();
  vnl_matrix<double> rel = RelativeTransform(I, I);
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
      EXPECT_NEAR(rel(i, j), (i == j ? 1.0 : 0.0), kEps);
}

// T_src = I, T_tgt = trans(1,2,3)  →  T_rel = T_tgt^{-1} = trans(-1,-2,-3)
TEST(RelativeTransform, SrcIdentityTgtTranslation) {
  vnl_matrix<double> rel = RelativeTransform(Identity44(), Trans(1.0, 2.0, 3.0));
  EXPECT_NEAR(rel(0, 3), -1.0, kEps);
  EXPECT_NEAR(rel(1, 3), -2.0, kEps);
  EXPECT_NEAR(rel(2, 3), -3.0, kEps);
  // Rotation block should be identity.
  for (int i = 0; i < 3; ++i)
    EXPECT_NEAR(rel(i, i), 1.0, kEps);
}

// T_rel = T_tgt^{-1} * T_src  →  T_tgt * T_rel should equal T_src exactly.
TEST(RelativeTransform, InverseProperty) {
  vnl_matrix<double> T_src = RotZ(30.0);
  T_src(0, 3) = 1.0;  T_src(1, 3) = 2.0;  T_src(2, 3) = 0.5;
  vnl_matrix<double> T_tgt = RotX(45.0);
  T_tgt(0, 3) = 0.5;  T_tgt(1, 3) = -1.0;  T_tgt(2, 3) = 3.0;

  vnl_matrix<double> T_rel = RelativeTransform(T_src, T_tgt);
  vnl_matrix<double> recovered = T_tgt * T_rel;
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
      EXPECT_NEAR(recovered(i, j), T_src(i, j), 1e-10)
          << "at (" << i << "," << j << ")";
}

// ── RotationErrorDeg ──────────────────────────────────────────────────────────

TEST(RotationErrorDeg, IdenticalTransforms) {
  vnl_matrix<double> T = RotZ(30.0);
  EXPECT_NEAR(RotationErrorDeg(T, T), 0.0, kEps);
}

TEST(RotationErrorDeg, Known90Degrees) {
  EXPECT_NEAR(RotationErrorDeg(RotZ(90.0), Identity44()), 90.0, 1e-10);
}

TEST(RotationErrorDeg, Known180Degrees) {
  EXPECT_NEAR(RotationErrorDeg(RotX(180.0), Identity44()), 180.0, 1e-10);
}

// Error is symmetric: RotationErrorDeg(A, B) == RotationErrorDeg(B, A).
TEST(RotationErrorDeg, Symmetry) {
  vnl_matrix<double> T1 = RotZ(30.0);
  vnl_matrix<double> T2 = RotX(45.0);
  EXPECT_NEAR(RotationErrorDeg(T1, T2), RotationErrorDeg(T2, T1), kEps);
}

// Translation columns must not affect the rotation error.
TEST(RotationErrorDeg, TranslationIgnored) {
  vnl_matrix<double> T1 = RotZ(90.0);
  vnl_matrix<double> T2 = Identity44();
  T1(0, 3) = 5.0;  T1(1, 3) = -3.0;
  T2(0, 3) = 2.0;  T2(2, 3) =  7.0;
  EXPECT_NEAR(RotationErrorDeg(T1, T2), 90.0, 1e-10);
}

// ── RotationMagnitudeDeg ──────────────────────────────────────────────────────

TEST(RotationMagnitudeDeg, Identity) {
  EXPECT_NEAR(RotationMagnitudeDeg(Identity44()), 0.0, kEps);
}

TEST(RotationMagnitudeDeg, NinetyDegAroundZ) {
  EXPECT_NEAR(RotationMagnitudeDeg(RotZ(90.0)), 90.0, 1e-10);
}

TEST(RotationMagnitudeDeg, OneEightyDegAroundX) {
  EXPECT_NEAR(RotationMagnitudeDeg(RotX(180.0)), 180.0, 1e-10);
}

TEST(RotationMagnitudeDeg, FortyFiveDeg) {
  EXPECT_NEAR(RotationMagnitudeDeg(RotZ(45.0)), 45.0, 1e-10);
}

// TranslationIgnored: adding a translation to identity should still give 0°.
TEST(RotationMagnitudeDeg, TranslationIgnored) {
  vnl_matrix<double> T = Trans(10.0, -5.0, 3.0);
  EXPECT_NEAR(RotationMagnitudeDeg(T), 0.0, kEps);
}
