#include "gmmreg_base.h"

#include <chrono>

#include <assert.h>
#include <fstream>
#include <iostream>

#include "utils/io_utils.h"
#include "utils/match_utils.h"
#include "utils/misc_utils.h"
#include "utils/normalization_utils.h"

namespace gmmreg {

typedef std::chrono::high_resolution_clock::time_point TimeVar;

#define duration(a) \
  std::chrono::duration_cast<std::chrono::microseconds>(a).count()
#define time_now() std::chrono::high_resolution_clock::now()


void Base::Run(const char* f_config) {
  Initialize(f_config);
  vnl_vector<double> params;
  TimeVar t1 = time_now();
  StartRegistration(params);
  elapsed_time_in_ms_ = duration(time_now() - t1) / 1000.0;
  std::cout << "Registration took " << elapsed_time_in_ms_ << " milliseconds."
      << std::endl;
  SaveResults(f_config, params);
}

int Base::Initialize(const char* f_config) {
  if (PrepareInput(f_config) < 0) {
    return -1;
  }
  SetInitParams(f_config);
  PrepareCommonOptions(f_config);
  PrepareOwnOptions(f_config);
  TimeVar t1 = time_now();
  PrepareBasisKernel();
  initialization_time_in_ms_ = duration(time_now() - t1) / 1000.0;
  std::cout << "Init took " << initialization_time_in_ms_ << " milliseconds."
      << std::endl;
  return 0;
}

int Base::PrepareInput(const char* f_config) {
  char f_model[256] = {0}, f_scene[256] = {0};
  GetPrivateProfileString(common_section_, "model", NULL,
      f_model, 256, f_config);
  if (LoadMatrixFromTxt(f_model, model_) < 0) {
    return -1;
  }

  GetPrivateProfileString(common_section_, "scene", NULL,
      f_scene, 256, f_config);
  if (LoadMatrixFromTxt(f_scene, scene_) < 0) {
    return -1;
  }

  m_ = model_.rows();
  d_ = model_.cols();
  transformed_model_.set_size(m_, d_);
  s_ = scene_.rows();
  assert(scene_.cols() == d_);

  return 0;
}

int Base::SetCtrlPts(const char* filename) {
  if (strlen(filename) == 0) {
    std::cout << "The control point set is not specified, "
                 "the model points are used as control points." << std::endl;
    ctrl_pts_ = model_;
    n_ = ctrl_pts_.rows();
    return n_;
  } else {
    if (LoadMatrixFromTxt(filename, ctrl_pts_) < 0) {
      return -1;
    }
    assert(ctrl_pts_.cols() == d_);
    n_ = ctrl_pts_.rows();
    return n_;
  }
}

void Base::DenormalizeAll() {
  if (b_normalize_) {
    Denormalize(transformed_model_, scene_centroid_, scene_scale_);
    Denormalize(model_, scene_centroid_, scene_scale_);
    Denormalize(scene_, scene_centroid_, scene_scale_);
  }
}

void Base::SaveElaspedTime(const char* f_config) {
  vnl_vector<double> elapsed_time;
  elapsed_time.set_size(1);
  elapsed_time[0] = elapsed_time_in_ms_;
  char f_elasped_time[256] = {0};
  GetPrivateProfileString(this->common_section_, "elasped_time_in_ms", NULL,
                          f_elasped_time, 256, f_config);

  SaveVectorToAsciiFile(f_elasped_time, elapsed_time);
}

void Base::SaveTransformed(const char* filename,
    const vnl_vector<double>& params, const char* f_config) {
  std::ofstream outfile(filename, std::ios_base::out);
  PerformTransform(params);
  DenormalizeAll();
  transformed_model_.print(outfile);

  char section_correspondence[256] = "CORRESPONDENCE";
  int num = GetPrivateProfileInt(section_correspondence,
                                 "num_of_thresholds", 0, f_config);
  if (num > 0) {
    char s_min[256], s_max[256], s_pairs[256];
    GetPrivateProfileString(section_correspondence,
        "min_threshold", NULL, s_min, 255, f_config);
    GetPrivateProfileString(section_correspondence,
        "max_threshold", NULL, s_max, 255, f_config);
    GetPrivateProfileString(section_correspondence,
        "matched_pairs", NULL, s_pairs, 255, f_config);
    std::ofstream f_pair(s_pairs, std::ios_base::out);
    double min_threshold, max_threshold, interval;
    min_threshold = atof(s_min);
    max_threshold = atof(s_max);
    if (num == 1) {
      interval = 0.0f;
    } else {
      interval = (max_threshold - min_threshold) / (num - 1);
    }
    //vnl_matrix<double> working_M, working_S;
    vnl_matrix<double> dist;
    vnl_matrix<int> pairs;
    ComputeSquaredDistanceMatrix(transformed_model_, scene_, dist);
    for (int i = 0; i < num; ++i) {
      double threshold  = min_threshold + i * interval;
      // int n_match = find_working_pair(model, scene, transformed_model,
      //                                threshold, working_M, working_S);
      FindNearestPairs<double>(dist, pairs, threshold * threshold);
      //printf("%f : %d\n",threshold, n_match);
      f_pair << "distance threshold : " << threshold << std::endl;
      f_pair << "# of matched point pairs : " << pairs.cols() << std::endl;
      int j;
      for (j = 0; j < pairs.cols(); ++j) {
        f_pair.width(6);
        f_pair << std::left << pairs(0, j);
      }
      f_pair << std::endl;
      for (j = 0; j < pairs.cols(); ++j) {
        f_pair.width(6);
        f_pair << std::left << pairs(1, j);
      }
      f_pair << std::endl;
    }
  }
  std::cout << "Please find the transformed model set in "
            << filename << std::endl;
}

void Base::PrepareCommonOptions(const char* f_config) {
  b_normalize_ = GetPrivateProfileInt(section_, "normalize", 1, f_config);
  if (b_normalize_) {
    Normalize(model_, model_centroid_, model_scale_);
    Normalize(scene_, scene_centroid_, scene_scale_);
    Normalize(ctrl_pts_, model_centroid_, model_scale_);
  }
#ifdef USE_KDTREE
  scene_tree_.reset(new NanoflannTree<double>(scene_));
  scene_tree_->tree.buildIndex();
  model_tree_.reset(new NanoflannTree<double>(model_));
  model_tree_->tree.buildIndex();
#endif
}

void Base::MultiScaleOptions(const char* f_config) {
  level_ = GetPrivateProfileInt(section_, "level", 1, f_config);
  char s_scale[256] = {0}, s_func_evals[256] = {0};
  char delims[] = " -,;";
  GetPrivateProfileString(section_, "sigma", NULL, s_scale, 255, f_config);
  utils::parse_tokens(s_scale, delims, v_scale_);
  if (v_scale_.size() < level_) {
    std::cerr << " Need more 'sigma' parameters. " << std::endl;
    exit(1);
  }
  GetPrivateProfileString(section_, "max_function_evals", NULL,
      s_func_evals, 255, f_config);
  utils::parse_tokens(s_func_evals, delims, v_func_evals_);
  if (v_func_evals_.size() < level_) {
    std::cerr << " Need more 'max_func_evals' parameters. " << std::endl;
    exit(1);
  }
}

}  // namespace gmmreg
