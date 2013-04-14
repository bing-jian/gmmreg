#include <assert.h>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>

#include <vcl_iostream.h>
#include <vnl/algo/vnl_lbfgs.h>
#include <vnl/algo/vnl_qr.h>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_trace.h>

#include "gmmreg_utils.h"
#include "gmmreg_tps.h"

int gmmreg_tps::prepare_input(const char* f_config) {
  gmmreg_base::prepare_input(f_config);
  char f_ctrl_pts[256] = {0};
  GetPrivateProfileString(common_section, "ctrl_pts", NULL,
      f_ctrl_pts, 256, f_config);
  if (set_ctrl_pts(f_ctrl_pts) < 0) {
    //todo: compute the ctrl pts on the fly
    return -1;
  }
  return 0;
}

void gmmreg_tps::start_registration(vnl_vector<double>& params) {
  vnl_lbfgs minimizer(*func);
  func->set_gmmreg(this);
  for (unsigned int k = 0; k < level; ++k) {
    func->set_scale(v_scale[k]);
    func->set_lambda(v_lambda[k]);
    bool b_fix_affine = (v_affine[k] == 1);
    func->set_fix_affine(b_fix_affine);
    func->prepare_param_gradient();
    set_param(params);
    int n_max_func_evals = v_func_evals[k];
    minimizer.set_max_function_evals(n_max_func_evals);
    // For more options, see
    // http://public.kitware.com/vxl/doc/release/core/vnl/html/vnl__nonlinear__minimizer_8h-source.html
    minimizer.minimize(params);
    if (minimizer.get_failure_code() < 0) {
      break;
    }
  }
  vcl_cout << "registration done" << vcl_endl;
}

int gmmreg_tps::set_init_params(const char* f_config) {
  char f_init_affine[256] = {0}, f_init_tps[256] = {0};
  GetPrivateProfileString(common_section, "init_affine", NULL,
      f_init_affine, 256, f_config);
  set_init_affine(f_init_affine);
  GetPrivateProfileString(common_section, "init_tps", NULL,
      f_init_tps, 256, f_config);
  set_init_tps(f_init_tps);
  param_all.set_size(n, d);
  param_all.update(param_affine);
  param_all.update(param_tps, d + 1);
  return 0;
}

int gmmreg_tps::set_init_affine(const char* filename) {
  if (strlen(filename) == 0) {
    // set default affine parameters from identity transform
    assert(d > 0);
    param_affine.set_size(d + 1, d);
    // the first row is for translation
    param_affine.fill(0);
    // the next dxd matrix is for affine matrix
    vnl_matrix<double> id;
    id.set_size(d, d);
    id.set_identity();
    param_affine.update(id, 1);
    return 0;
  } else {
    std::ifstream infile(filename);
    param_affine.read_ascii(infile);
    assert(param_affine.cols() == d);
    assert(param_affine.rows() == (d + 1));
    return 1;
  }
}

int gmmreg_tps::set_init_tps(const char* filename) {
  if (strlen(filename) == 0) {
    assert(n - d - 1 > 0);
    assert(d > 0);
    param_tps.set_size(n - d - 1, d);
    param_tps.fill(0);
    return 0;
  } else {
    std::ifstream infile(filename, std::ios_base::in);
    param_tps.read_ascii(infile);
    assert(param_tps.cols() == d);
    assert(param_tps.rows() == (n - d - 1));
    return 1;
  }
}

void gmmreg_tps::prepare_basis_kernel() {
  //todo: detect singularity of the data
  vnl_matrix<double> K, U;
  ComputeTPSKernel(model, ctrl_pts, U, K);
  m = model.rows();
  vnl_matrix<double> Pm;
  Pm.set_size(m, d + 1);
  Pm.set_column(0, 1);
  Pm.update(model, 0, 1);

  vnl_matrix<double> Pn;
  Pn.set_size(n, d + 1);
  Pn.set_column(0, 1);
  Pn.update(ctrl_pts, 0, 1);
  /* should use SVD(Pn), but vnl's SVD is an ``economy-size'' SVD  */

  /* vnl_svd<double> SVD(Pn.transpose());
  vnl_matrix<double> VV = SVD.V();
  std::cout << VV.rows() << " " << VV.cols() << std::endl;
  save_matrix("./VV.txt", VV);
  */

  vnl_qr<double> qr(Pn);
  vnl_matrix<double> V = qr.Q();
  vnl_matrix<double> PP = V.extract(n, n - d - 1, 0, d + 1);
  basis.set_size(m, n);
  basis.update(Pm);
  basis.update(U * PP, 0, d + 1);
  kernel = PP.transpose() * K * PP;
}

void gmmreg_tps::perform_transform(const vnl_vector<double> &x) {
  set_affine_and_tps(x);
  transformed_model = basis * param_all;
  save_matrix("./param_all.txt", param_all);
  save_matrix("./basis.txt", basis);
}

double gmmreg_tps::bending_energy() {
  return vnl_trace(param_tps.transpose() * kernel * param_tps);
}

void gmmreg_tps::compute_gradient(double lambda,
    const vnl_matrix<double>& gradient, vnl_matrix<double>& grad_all) {
  grad_all.fill(0);
  if (lambda > 0) {
    grad_all.update(2 * lambda * kernel * param_tps, d+1);
  }
  grad_all += basis.transpose() * gradient;
}

void gmmreg_tps::save_results(const char* f_config,
    const vnl_vector<double>& params) {
  char f_transformed[256] = {0};
  char f_final_affine[256] = {0};
  char f_final_tps[256] = {0};

  GetPrivateProfileString(common_section, "final_affine", NULL,
      f_final_affine, 255, f_config);
  GetPrivateProfileString(common_section, "final_tps", NULL,
      f_final_tps, 255, f_config);
  GetPrivateProfileString(common_section, "transformed_model", NULL,
      f_transformed, 255, f_config);

  save_transformed(f_transformed, params, f_config);
  save_matrix(f_final_affine, param_affine);
  save_matrix(f_final_tps, param_tps);
}


void gmmreg_tps::prepare_own_options(const char* f_config) {
  multi_scale_options(f_config);
  char delims[] = " -,;";
  char s_lambda[256] = {0};
  GetPrivateProfileString(section, "lambda", NULL, s_lambda, 255, f_config);
  parse_tokens(s_lambda, delims, v_lambda);
  if (v_lambda.size() < level) {
    std::cerr<< " too many levels " << std::endl;
    exit(1);
  }
  char s_affine[256] = {0};
  GetPrivateProfileString(section, "fix_affine", NULL,
      s_affine, 255, f_config);
  parse_tokens(s_affine, delims, v_affine);
  if (v_affine.size() < level) {
    std::cerr<< " too many levels " << std::endl;
    exit(1);
  }
}

void gmmreg_tps::set_param(vnl_vector<double>& x0) {
  int k = 0;
  x0.set_size(func->get_number_of_unknowns());
  x0.fill(0);
  if (!func->fix_affine) { // x0 includes affine
    for (int i = 0; i < param_affine.rows(); ++i) {
      for (int j = 0; j < d; ++j, ++k) {
        x0[k] = param_affine(i, j);
      }
    }
  }
  for (int i = 0; i < param_tps.rows(); ++i) {
    for (int j = 0; j < d; ++j, ++k) {
      x0[k] = param_tps(i, j);
    }
  }
}

void gmmreg_tps::set_affine_and_tps(const vnl_vector<double>& x) {
  /* reshape x, assuming x is row major; */
  int rows_x = x.size() / d;
  if (func->fix_affine) {   // affine is given, x does not include affine
    param_all.update(param_affine);
    for (int i = 0, k = 0; i < rows_x; ++i) {
      for (int j = 0; j < d; ++j, ++k) {
        param_tps(i, j) = x[k];
      }
    }
    param_all.update(param_tps, d + 1);
  } else { // affine is not given, x includes affine already
    for (int i = 0, k = 0; i < rows_x; ++i) {
      for (int j = 0; j < d; ++j, ++k) {
        param_all(i, j) = x[k];
      }
    }
    param_affine = param_all.extract(d + 1, d);
    param_tps = param_all.extract(rows_x - d - 1, d, d + 1);
  }
}
