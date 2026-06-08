// Sequential pairwise 3-D rigid registration over a directory of point clouds.
//
// Loads consecutively numbered point clouds (pcloud_XXXXX.txt), registers
// each pair (src, src+step) using gmmreg::RigidRegistration, and writes the
// recovered 4×4 transform (and, optionally, the transformed point cloud) to
// output_dir.  A summary CSV is always written to output_dir/summary.csv.
//
// Optimization parameters (normalize, sigma, max_function_evals) are read
// from a config file's [GMMREG_OPT] section — the same INI format used by
// gmmreg_demo.  See lounge.ini for an example.
//
// When --gt_dir is given the program also loads the corresponding GT poses
// (gt_pose_XXXXX.txt, 4×4 absolute poses) and reports the geodesic rotation
// error between the estimated relative transform and the GT relative transform.
//
// Usage
// -----
//   seq_rigid3d \
//       --pcloud_dir <dir>             \
//       --output_dir <dir>             \
//       --config     <ini_file>        \
//     [ --gt_dir <dir>               ] \
//     [ --step           <int>       ] \   (default: 1)
//     [ --start          <int>       ] \   (1-based frame index, default: 1)
//     [ --end            <int>       ] \   (default: auto-detected from dir)
//     [ --save_transformed <0|1>     ]     (default: 0)

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

#include <vnl/vnl_matrix.h>

#ifdef WIN32
#include <windows.h>
#else
#include "port_ini.h"
#endif

#include "gmmreg_rigid.h"
#include "utils/io_utils.h"
#include "utils/misc_utils.h"
#include "utils/transform_utils.h"

namespace fs = std::filesystem;

// ── Config ────────────────────────────────────────────────────────────────────

struct Config {
  std::string pcloud_dir;
  std::string output_dir;
  std::string config_file;
  std::string gt_dir;
  int  step             = 1;
  int  start            = 1;   // 1-based
  int  end              = -1;  // -1 = auto-detect
  int  max_pairs        = 10;   // -1 = no limit
  bool save_transformed = false;
  bool dry_run          = false;
};

struct RegConfig {
  bool normalize = true;
  std::vector<double> scales        = {1.0};
  std::vector<int>    max_func_evals = {100};
};

// Read optimization parameters from [GMMREG_OPT] section of an INI file.
// Mirrors Base::MultiScaleOptions: parse all tokens then truncate to level.
static RegConfig LoadRegConfig(const char* ini_file) {
  RegConfig rc;
  static const char* kSection = "GMMREG_OPT";
  char delims[] = " -,;";

  rc.normalize = (GetPrivateProfileInt(kSection, "normalize", 1, ini_file) != 0);

  int level = GetPrivateProfileInt(kSection, "level", 1, ini_file);

  char s_scale[256] = {0};
  GetPrivateProfileString(kSection, "sigma", "1.0", s_scale, 255, ini_file);
  std::vector<float> v_scale;
  gmmreg::utils::parse_tokens(s_scale, delims, v_scale);
  if (static_cast<int>(v_scale.size()) < level) {
    std::cerr << "ERROR: need at least " << level
              << " sigma values in [GMMREG_OPT]\n";
    std::exit(1);
  }
  v_scale.resize(level);
  rc.scales.assign(v_scale.begin(), v_scale.end());

  char s_evals[256] = {0};
  GetPrivateProfileString(kSection, "max_function_evals", "100", s_evals, 255,
                          ini_file);
  std::vector<int> v_evals;
  gmmreg::utils::parse_tokens(s_evals, delims, v_evals);
  if (static_cast<int>(v_evals.size()) < level) {
    std::cerr << "ERROR: need at least " << level
              << " max_function_evals values in [GMMREG_OPT]\n";
    std::exit(1);
  }
  v_evals.resize(level);
  rc.max_func_evals = v_evals;

  return rc;
}

// ── Path helpers ──────────────────────────────────────────────────────────────

static std::string PcloudPath(const std::string& dir, int n) {
  std::ostringstream s;
  s << dir << "/pcloud_" << std::setw(5) << std::setfill('0') << n << ".txt";
  return s.str();
}

static std::string GtPosePath(const std::string& dir, int n) {
  std::ostringstream s;
  s << dir << "/gt_pose_" << std::setw(5) << std::setfill('0') << n << ".txt";
  return s.str();
}

static std::string PairTag(int src, int tgt) {
  std::ostringstream s;
  s << std::setw(5) << std::setfill('0') << src
    << "_" << std::setw(5) << std::setfill('0') << tgt;
  return s.str();
}

static int DetectLastFrame(const std::string& dir) {
  int last = 0;
  for (const auto& entry : fs::directory_iterator(dir)) {
    if (!entry.is_regular_file()) continue;
    const std::string name = entry.path().filename().string();
    if (name.rfind("pcloud_", 0) != 0) continue;
    try {
      int n = std::stoi(name.substr(7, 5));
      last = std::max(last, n);
    } catch (...) {}
  }
  return last;
}

// ── Stats helpers ─────────────────────────────────────────────────────────────

struct Stats {
  double avg, min, max, median;
};

static Stats ComputeStats(std::vector<double> v) {
  Stats s{};
  if (v.empty()) return s;
  std::sort(v.begin(), v.end());
  s.min    = v.front();
  s.max    = v.back();
  s.avg    = std::accumulate(v.begin(), v.end(), 0.0) / v.size();
  size_t n = v.size();
  s.median = (n % 2 == 0) ? (v[n / 2 - 1] + v[n / 2]) * 0.5 : v[n / 2];
  return s;
}

// Write a simple JSON results file to output_path.
static void WriteJsonReport(const std::string& output_path,
                            int n_ok, int n_failed,
                            int start, int end, int step,
                            const std::vector<double>& times_ms,
                            const std::vector<double>& rot_errors,
                            const std::vector<double>& gt_diffs) {
  std::ofstream f(output_path);
  auto fmt4 = [](double v) {
    std::ostringstream s;
    s << std::fixed << std::setprecision(4) << v;
    return s.str();
  };

  f << "{\n"
    << "  \"pairs_ok\": "     << n_ok     << ",\n"
    << "  \"pairs_failed\": " << n_failed << ",\n"
    << "  \"frames\": {\"start\": " << start
    <<               ", \"end\": "   << end
    <<               ", \"step\": "  << step << "},\n";

  if (!times_ms.empty()) {
    Stats t = ComputeStats(times_ms);
    f << "  \"runtime_ms\": {"
      << "\"avg\": " << fmt4(t.avg) << ", "
      << "\"min\": " << fmt4(t.min) << ", "
      << "\"max\": " << fmt4(t.max) << ", "
      << "\"median\": " << fmt4(t.median) << "},\n"
      << "  \"fps\": {"
      << "\"avg\": "    << fmt4(1000.0 / t.avg)    << ", "
      << "\"median\": " << fmt4(1000.0 / t.median) << "},\n";
  }

  if (!rot_errors.empty()) {
    Stats r = ComputeStats(rot_errors);
    f << "  \"rot_error_deg\": {"
      << "\"avg\": " << fmt4(r.avg) << ", "
      << "\"min\": " << fmt4(r.min) << ", "
      << "\"max\": " << fmt4(r.max) << ", "
      << "\"median\": " << fmt4(r.median) << "}";
    if (!gt_diffs.empty()) f << ",";
    f << "\n";
  }

  if (!gt_diffs.empty()) {
    Stats g = ComputeStats(gt_diffs);
    f << "  \"gt_pose_diff_deg\": {"
      << "\"avg\": " << fmt4(g.avg) << ", "
      << "\"min\": " << fmt4(g.min) << ", "
      << "\"max\": " << fmt4(g.max) << ", "
      << "\"median\": " << fmt4(g.median) << "}\n";
  }

  f << "}\n";
}

// ── Entry point ───────────────────────────────────────────────────────────────

int main(int argc, char* argv[]) {
  Config cfg;

  for (int i = 1; i + 1 < argc; i += 2) {
    std::string key = argv[i];
    std::string val = argv[i + 1];
    if      (key == "--pcloud_dir")       cfg.pcloud_dir       = val;
    else if (key == "--output_dir")       cfg.output_dir       = val;
    else if (key == "--config")           cfg.config_file      = val;
    else if (key == "--gt_dir")           cfg.gt_dir           = val;
    else if (key == "--step")             cfg.step             = std::stoi(val);
    else if (key == "--start")            cfg.start            = std::stoi(val);
    else if (key == "--end")              cfg.end              = std::stoi(val);
    else if (key == "--max_pairs")        cfg.max_pairs        = std::stoi(val);
    else if (key == "--save_transformed") cfg.save_transformed = (val != "0");
    else if (key == "--dry_run")          cfg.dry_run          = (val != "0");
    else { std::cerr << "Unknown option: " << key << "\n"; return 1; }
  }

  // --dry_run only requires --pcloud_dir; skip --output_dir and --config check.
  if (cfg.pcloud_dir.empty() ||
      (!cfg.dry_run && (cfg.output_dir.empty() || cfg.config_file.empty()))) {
    std::cerr
        << "Usage: seq_rigid3d\n"
        << "  --pcloud_dir <dir>            point cloud input directory\n"
        << "  --output_dir <dir>            output directory\n"
        << "  --config     <ini_file>       INI file with [GMMREG_OPT] section\n"
        << "  [--gt_dir <dir>]              GT pose directory (for error)\n"
        << "  [--step <n>]                  frame step size (default: 1)\n"
        << "  [--start <n>]                 first frame, 1-based (default: 1)\n"
        << "  [--end <n>]                   last frame (default: auto)\n"
        << "  [--max_pairs <n>]             stop after n pairs; -1 = no limit (default: 10)\n"
        << "  [--save_transformed <0|1>]    save transformed point cloud (default: 0)\n"
        << "  [--dry_run <0|1>]             check data availability only (default: 0)\n";
    return 1;
  }

  if (cfg.end < 0)
    cfg.end = DetectLastFrame(cfg.pcloud_dir);

  if (cfg.end <= 0) {
    std::cerr << "ERROR: no pcloud_XXXXX.txt files found in " << cfg.pcloud_dir << "\n";
    return 1;
  }

  const bool have_gt = !cfg.gt_dir.empty();

  if (have_gt && !fs::is_directory(cfg.gt_dir)) {
    std::cerr << "ERROR: gt_dir does not exist: " << cfg.gt_dir << "\n";
    return 1;
  }

  // ── Dry-run: check data availability + config, report, then exit ───────────
  if (cfg.dry_run) {
    // Config file (optional in dry-run, but validated and printed if given).
    bool config_ok = false;
    if (!cfg.config_file.empty()) {
      if (!fs::exists(cfg.config_file)) {
        std::cerr << "ERROR: config file not found: " << cfg.config_file << "\n";
      } else {
        const RegConfig rc = LoadRegConfig(cfg.config_file.c_str());
        config_ok = true;
        std::cout << "Config: " << cfg.config_file << "\n"
                  << "  normalize      : " << rc.normalize << "\n"
                  << "  levels         : " << rc.scales.size() << "\n"
                  << "  sigma          :";
        for (double s : rc.scales) std::cout << " " << s;
        std::cout << "\n  max_func_evals :";
        for (int e : rc.max_func_evals) std::cout << " " << e;
        std::cout << "\n";
      }
    } else {
      std::cout << "Config: (none provided)\n";
    }

    // Data availability scan.
    int n_pairs = 0, n_missing_pcloud = 0, n_missing_gt = 0;
    for (int src = cfg.start; src + cfg.step <= cfg.end; ++src) {
      if (cfg.max_pairs >= 0 && n_pairs >= cfg.max_pairs) break;
      const int tgt = src + cfg.step;
      ++n_pairs;
      if (!fs::exists(PcloudPath(cfg.pcloud_dir, src)) ||
          !fs::exists(PcloudPath(cfg.pcloud_dir, tgt)))
        ++n_missing_pcloud;
      if (have_gt &&
          (!fs::exists(GtPosePath(cfg.gt_dir, src)) ||
           !fs::exists(GtPosePath(cfg.gt_dir, tgt))))
        ++n_missing_gt;
    }
    std::cout << "Data summary\n"
              << "  pcloud_dir     : " << cfg.pcloud_dir << "\n"
              << "  gt_dir         : " << (have_gt ? cfg.gt_dir : "(none)") << "\n"
              << "  frames         : " << cfg.start << " .. " << cfg.end
              << "  step=" << cfg.step << "\n"
              << "  max_pairs      : " << (cfg.max_pairs < 0 ? "unlimited" : std::to_string(cfg.max_pairs)) << "\n"
              << "  pairs to check : " << n_pairs << "\n"
              << "  pcloud missing : " << n_missing_pcloud << "\n";
    if (have_gt)
      std::cout << "  gt_pose missing: " << n_missing_gt << "\n";

    const bool data_ok = (n_missing_pcloud == 0 && n_missing_gt == 0);
    std::cout << (data_ok && (!cfg.config_file.empty() ? config_ok : true)
                      ? "Ready to run.\n"
                      : "NOT ready — fix errors above before running.\n");
    return (data_ok && (!cfg.config_file.empty() ? config_ok : true)) ? 0 : 1;
  }
  // ───────────────────────────────────────────────────────────────────────────

  // Non-dry-run: --config and --output_dir are required.
  const RegConfig rc = LoadRegConfig(cfg.config_file.c_str());

  // Early stop if a needed GT file is absent (spot-check first and last).
  if (have_gt) {
    for (int n : {cfg.start, cfg.end}) {
      std::string p = GtPosePath(cfg.gt_dir, n);
      if (!fs::exists(p)) {
        std::cerr << "ERROR: GT pose file not found: " << p << "\n";
        return 1;
      }
    }
  }

  fs::create_directories(cfg.output_dir);

  std::ofstream summary(cfg.output_dir + "/summary.csv");
  summary << "src,tgt,elapsed_ms";
  if (have_gt) summary << ",rot_error_deg";
  summary << "\n";

  int n_ok = 0, n_failed = 0;
  std::vector<double> v_times_ms, v_rot_errors, v_gt_diffs;

  for (int src = cfg.start; src + cfg.step <= cfg.end; ++src) {
    if (cfg.max_pairs >= 0 && (n_ok + n_failed) >= cfg.max_pairs) break;
    const int tgt = src + cfg.step;

    vnl_matrix<double> model, scene;
    if (gmmreg::LoadMatrixFromTxt(PcloudPath(cfg.pcloud_dir, src).c_str(), model) < 0 ||
        gmmreg::LoadMatrixFromTxt(PcloudPath(cfg.pcloud_dir, tgt).c_str(), scene) < 0) {
      std::cerr << "  WARNING: cannot load pair " << src << "/" << tgt << "\n";
      ++n_failed;
      continue;
    }

    gmmreg::RegistrationInput input;
    input.model          = model;
    input.scene          = scene;
    input.normalize      = rc.normalize;
    input.scales         = rc.scales;
    input.max_func_evals = rc.max_func_evals;

    gmmreg::RigidRegistration reg;
    const auto t0 = std::chrono::high_resolution_clock::now();
    if (reg.RunWithData(input) < 0) {
      std::cerr << "  WARNING: RunWithData failed for pair " << src << "/" << tgt << "\n";
      ++n_failed;
      continue;
    }
    const double elapsed_ms = std::chrono::duration<double, std::milli>(
        std::chrono::high_resolution_clock::now() - t0).count();
    v_times_ms.push_back(elapsed_ms);

    const std::string tag = PairTag(src, tgt);

    vnl_matrix<double> matrix;
    reg.GetTransformMatrix(&matrix);
    gmmreg::SaveMatrixToAsciiFile(
        (cfg.output_dir + "/transform_" + tag + ".txt").c_str(), matrix);

    if (cfg.save_transformed)
      gmmreg::SaveMatrixToAsciiFile(
          (cfg.output_dir + "/transformed_" + tag + ".txt").c_str(),
          reg.GetTransformedModel());

    summary << src << "," << tgt
            << "," << std::fixed << std::setprecision(3) << elapsed_ms;

    if (have_gt) {
      vnl_matrix<double> gt_src, gt_tgt;
      if (gmmreg::LoadMatrixFromTxt(GtPosePath(cfg.gt_dir, src).c_str(), gt_src) > 0 &&
          gmmreg::LoadMatrixFromTxt(GtPosePath(cfg.gt_dir, tgt).c_str(), gt_tgt) > 0) {
        const vnl_matrix<double> T_gt_rel = gmmreg::RelativeTransform(gt_src, gt_tgt);
        const double rot_err = gmmreg::RotationErrorDeg(matrix, T_gt_rel);
        v_rot_errors.push_back(rot_err);
        v_gt_diffs.push_back(gmmreg::RotationMagnitudeDeg(T_gt_rel));
        summary << "," << std::setprecision(4) << rot_err;
      } else {
        summary << ",nan";
      }
    }
    summary << "\n";

    ++n_ok;
    if (n_ok % 100 == 0)
      std::cout << "  " << n_ok << " pairs done\n";
  }

  const std::string json_path = cfg.output_dir + "/results.json";
  WriteJsonReport(json_path, n_ok, n_failed,
                  cfg.start, cfg.end, cfg.step,
                  v_times_ms, v_rot_errors, v_gt_diffs);

  std::cout << "Done.  ok=" << n_ok << "  failed=" << n_failed << "\n";
  std::cout << "Results in " << cfg.output_dir << "/\n";
  std::cout << "JSON report: " << json_path << "\n";
  return (n_failed > 0) ? 1 : 0;
}
