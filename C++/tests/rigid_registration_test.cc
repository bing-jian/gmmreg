#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

#include <gtest/gtest.h>
#include <vnl/vnl_matrix.h>

#include "gmmreg_api.h"
#include "gmmreg_rigid.h"
#include "utils/io_utils.h"
#include "utils/json_utils.h"

// ── Registration constants (easy to adjust) ───────────────────────────────────

constexpr int    kNormalize    = 1;
constexpr double kSigma        = 1.0;
constexpr int    kMaxFuncEvals = 100;
constexpr double kRmseThreshold = 0.5;  // original-coordinate units

// ── Helpers ───────────────────────────────────────────────────────────────────

namespace {

constexpr double kPi = 3.14159265358979323846;

// Apply 2D rigid: out[i] = pts[i] * R.T + [tx, ty].
// Matches RigidRegistration::PerformTransform with params=[tx, ty, theta].
void ApplyRigidTransform(const vnl_matrix<double>& pts,
                         double tx, double ty, double theta,
                         vnl_matrix<double>* out) {
  double c = std::cos(theta), s = std::sin(theta);
  int m = pts.rows();
  out->set_size(m, 2);
  for (int i = 0; i < m; ++i) {
    (*out)(i, 0) = pts(i, 0) * c - pts(i, 1) * s + tx;
    (*out)(i, 1) = pts(i, 0) * s + pts(i, 1) * c + ty;
  }
}

// Per-point RMSE between a flat row-major matrix and a vnl_matrix reference.
double ComputeRMSE(const std::vector<double>& flat,
                   const vnl_matrix<double>& ref) {
  int m = ref.rows(), d = ref.cols();
  double sum = 0.0;
  for (int i = 0; i < m; ++i)
    for (int j = 0; j < d; ++j) {
      double e = flat[static_cast<size_t>(i * d + j)] - ref(i, j);
      sum += e * e;
    }
  return std::sqrt(sum / m);
}

// Per-point RMSE between two vnl_matrices of the same shape.
double ComputeRMSE(const vnl_matrix<double>& a, const vnl_matrix<double>& b) {
  double sum = 0.0;
  for (int i = 0; i < a.rows(); ++i)
    for (int j = 0; j < a.cols(); ++j) {
      double e = a(i, j) - b(i, j);
      sum += e * e;
    }
  return std::sqrt(sum / a.rows());
}

// ── Parameterized test fixture ─────────────────────────────────────────────────

struct RigidParams {
  double tx, ty, theta;
  const char* name;
};

class RigidRegistrationTest
    : public ::testing::TestWithParam<RigidParams> {
 protected:
  void SetUp() override {
    const RigidParams& p = GetParam();
    prefix_       = std::string("rigid_") + p.name;
    f_model_      = prefix_ + "_model.txt";
    f_config_     = prefix_ + ".ini";
    f_transformed_= prefix_ + "_transformed.txt";
    f_output_     = prefix_ + "_output.json";
    f_final_      = prefix_ + "_final.txt";

    ASSERT_GT(gmmreg::LoadMatrixFromTxt("cmu_road_data.txt", scene_), 0)
        << "cmu_road_data.txt not found in working directory";

    ApplyRigidTransform(scene_, p.tx, p.ty, p.theta, &model_);
    gmmreg::SaveMatrixToAsciiFile(f_model_.c_str(), model_);
    WriteConfig();
  }

  void TearDown() override {
    for (const auto& f : {f_model_, f_config_, f_transformed_, f_output_, f_final_})
      std::remove(f.c_str());
  }

  void WriteConfig() {
    std::ofstream f(f_config_);
    f << "[FILES]\n"
      << "model = ./" << f_model_ << "\n"
      << "scene = ./cmu_road_data.txt\n"
      << "transformed_model = ./" << f_transformed_ << "\n"
      << "output_json = ./" << f_output_ << "\n"
      << "final_rigid = ./" << f_final_ << "\n"
      << "\n[GMMREG_OPT]\n"
      << "normalize = "          << kNormalize    << "\n"
      << "level = 1\n"
      << "sigma = "              << kSigma        << "\n"
      << "max_function_evals = " << kMaxFuncEvals << "\n";
  }

  std::string prefix_, f_model_, f_config_, f_transformed_, f_output_, f_final_;
  vnl_matrix<double> scene_, model_;
};

}  // namespace

// ── Test ──────────────────────────────────────────────────────────────────────

TEST_P(RigidRegistrationTest, RecoverRigidTransform) {
  ASSERT_EQ(gmmreg_api(f_config_.c_str(), "rigid"), 0)
      << "gmmreg_api failed for " << GetParam().name;

  std::string json = gmmreg::ReadAll(f_output_.c_str());
  ASSERT_FALSE(json.empty()) << "output JSON was not created";

  auto flat = gmmreg::ParseMatrix(json, "transformed_model");
  ASSERT_EQ(static_cast<int>(flat.size()), scene_.rows() * scene_.cols())
      << "transformed_model element count mismatch";

  double rmse = ComputeRMSE(flat, scene_);
  EXPECT_LT(rmse, kRmseThreshold)
      << "RMSE = " << rmse << " for transform '" << GetParam().name << "'";
}

INSTANTIATE_TEST_SUITE_P(
    CleanData2D,
    RigidRegistrationTest,
    ::testing::Values(
        RigidParams{ 0.5,  0.0,         0.0,       "translation_x"  },
        RigidParams{ 0.0,  0.5,         0.0,       "translation_y"  },
        RigidParams{ 0.0,  0.0,  kPi / 12.0,       "rotation_15deg" },
        RigidParams{ 0.3,  0.2,  kPi / 12.0,       "combined_small" }
    ),
    [](const ::testing::TestParamInfo<RigidParams>& info) {
      return std::string(info.param.name);
    }
);

// ── No-config tests: use RunWithData directly ─────────────────────────────────

class RigidRegistrationNoConfigTest
    : public ::testing::TestWithParam<RigidParams> {
 protected:
  void SetUp() override {
    ASSERT_GT(gmmreg::LoadMatrixFromTxt("cmu_road_data.txt", scene_), 0)
        << "cmu_road_data.txt not found in working directory";
    ApplyRigidTransform(scene_, GetParam().tx, GetParam().ty,
                        GetParam().theta, &model_);
  }
  vnl_matrix<double> scene_, model_;
};

TEST_P(RigidRegistrationNoConfigTest, RecoverRigidTransform) {
  gmmreg::RigidRegistration reg;
  ASSERT_EQ(reg.RunWithData(model_, scene_,
                             kNormalize != 0,
                             {kSigma},
                             {kMaxFuncEvals}), 0)
      << "RunWithData failed for " << GetParam().name;

  double rmse = ComputeRMSE(reg.GetTransformedModel(), scene_);
  EXPECT_LT(rmse, kRmseThreshold)
      << "RMSE = " << rmse << " for transform '" << GetParam().name << "'";
}

INSTANTIATE_TEST_SUITE_P(
    CleanData2D,
    RigidRegistrationNoConfigTest,
    ::testing::Values(
        RigidParams{ 0.5,  0.0,         0.0,       "translation_x"  },
        RigidParams{ 0.0,  0.5,         0.0,       "translation_y"  },
        RigidParams{ 0.0,  0.0,  kPi / 12.0,       "rotation_15deg" },
        RigidParams{ 0.3,  0.2,  kPi / 12.0,       "combined_small" }
    ),
    [](const ::testing::TestParamInfo<RigidParams>& info) {
      return std::string(info.param.name);
    }
);
