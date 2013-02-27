#include "gmmreg_cpd.h"

#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <vcl_iostream.h>
#include <vnl/algo/vnl_determinant.h>
#include <vnl/algo/vnl_qr.h>
#include <vnl/algo/vnl_svd.h>
#include <vnl/vnl_trace.h>

#include "gmmreg_utils.h"

int gmmreg_cpd::prepare_input(const char* f_config) {
  gmmreg_base::prepare_input(f_config);
  char f_ctrl_pts[80]={0};
  GetPrivateProfileString(common_section, "ctrl_pts", NULL,
      f_ctrl_pts, 80, f_config);
  if (set_ctrl_pts(f_ctrl_pts) < 0) {
    //todo: compute the ctrl pts on the fly
    return -1;
  }
  return 0;
}

void gmmreg_cpd_grbf::prepare_basis_kernel() {
  ComputeGaussianKernel(model, ctrl_pts, basis, kernel, beta);
  Gtranspose = basis.transpose();
  dPG.set_size(m, n);
  dPY0.set_size(m, d);
  //vnl_qr<double> qr(Gtranspose*dPG+lambda*sigma*sigma*kernel); //, 1e-18);
  //invG = qr.inverse()*Gtranspose;
}

double gmmreg_cpd_grbf::update_param() {
  double bending = vnl_trace(param_all.transpose()*kernel*param_all);
  double row_sum;
  for (int i = 0;i < m;++i) {
    row_sum  = P.get_row(i).sum();
    dPG.set_row(i, row_sum * basis.get_row(i));
    dPY0.set_row(i, row_sum * model.get_row(i));
  }
  vnl_qr<double> qr(Gtranspose*dPG+lambda*sigma*sigma*kernel); //, 1e-18);
  param_all = qr.solve(Gtranspose*(P*scene-dPY0));
  //param_all = invG*(P*scene-dPY0);
  return bending;
}

void gmmreg_cpd_tps::prepare_basis_kernel() {
  vnl_matrix<double> K, U;
  ComputeTPSKernel(model, ctrl_pts, U, K);

  //n = ctrl_pts.rows();
  vnl_matrix<double> Pn;
  Pn.set_size(n, d+1);
  Pn.set_column(0,1);
  Pn.update(ctrl_pts, 0, 1);
  vnl_qr<double> qr(Pn);
  vnl_matrix<double> V=qr.Q();
  vnl_matrix<double> PP = V.extract(n,n-d-1,0,d+1);
  kernel = PP.transpose()*K*PP;

  //m = model.rows();
  vnl_matrix<double> Pm;
  Pm.set_size(m, d+1);
  Pm.set_column(0,1);
  Pm.update(model, 0, 1);

  basis.set_size(m, n);
  basis.update(Pm);
  basis.update(U*PP,0,d+1);

  G = U*PP;
  vnl_matrix<double> Gtranspose;
  Gtranspose = G.transpose();

  vnl_qr<double> qr_Pm(Pm);
  Q1 = qr_Pm.Q().extract(m,d+1);
  Q2 = qr_Pm.Q().extract(m,m-d-1,0,d+1);
  R = qr_Pm.R().extract(d+1,d+1);

  vnl_svd<double> svd(R);
  invR = svd.inverse();
  nP.set_size(m,s);

  vnl_matrix<double> A = G.transpose()*Q2*Q2.transpose()*G+lambda*kernel;
  vnl_qr<double> qr2(A);
  invG = qr2.inverse()*G.transpose()*Q2*Q2.transpose();
}

double gmmreg_cpd_tps::update_param() {
  tps = param_all.extract(n-d-1,d,d+1,0);
  //std::cout << "before: tps " << tps.array_two_norm() << std::endl;
  double bending = vnl_trace(tps.transpose()*kernel*tps);
  //std::cout << "bending = " << bending << std::endl;
  double row_sum;
  for (int i = 0; i < m; ++i) {
    row_sum  = P.get_row(i).sum();
    if (row_sum > eps) {
      nP.set_row(i, P.get_row(i)/row_sum);
    }
  }
  //std::cout << "before: nP " << nP.array_two_norm() << std::endl;
  //vnl_qr<double> qr(G.transpose()*Q2*Q2.transpose()*G+lambda*kernel);
  //tps = qr.solve(G.transpose()*Q2*Q2.transpose()*(nP*scene));
  tps = invG*(nP*scene);
  affine = invR*Q1.transpose()*(nP*scene-model-G*tps);
  param_all.update(affine);
  param_all.update(tps,d+1);
  //std::cout << "after: tps" << tps.array_two_norm() << std::endl;
  return bending;
}

void gmmreg_cpd::start_registration(vnl_vector<double>& params) {
  int EMiter, iter = 0;
  double Eu, E_old;
  double E = 1;
  //outlier_term = outliers*pow((2*sigma*sigma*3.1415926),0.5*d);
  double ntol = tol + 10;
  vnl_matrix<double> dP;
  vnl_matrix<double> prev_model;
  vnl_matrix<double> moving_model(model);
  //vnl_matrix<double> eye(n,n);
  //eye.set_identity();
  /*prepare_basis_kernel(); already done */
  P.set_size(model.rows(), scene.rows());
  double bending;
  //column_sum.set_size(s);
  while ((iter < max_iter) && (ntol > tol)) {
    //std::cout << "iter=" << iter << "\n";
    EMiter = 0;   EMtol = tol + 10;
    prev_model = moving_model;
    while ((EMiter < max_em_iter) && (EMtol > tol)) {
      //std::cout << "EMiter="<<EMiter<< "\t";
      //std::cout << "E="<<E<<"\t";
      //std::cout << "sigma="<<sigma<<std::endl;
      E_old = E;
      compute_P(moving_model, scene, P, Eu, sigma, outliers);
      bending = update_param();
      moving_model = model + basis * param_all;
      E = Eu + (lambda / 2) * bending;
      EMtol = fabs(E_old - E) / E_old;
      ++EMiter;
    }
    sigma *= anneal;
    ++iter;
    ntol = (moving_model - prev_model).array_two_norm();
  }
}

int gmmreg_cpd::set_init_params(const char* f_config) {
  char f_init_params[80] = {0};
  GetPrivateProfileString("Files", "init_params", NULL,
      f_init_params, 80, f_config);
  if (strlen(f_init_params) == 0) {
    assert(n > 0);
    assert(d > 0);
    param_all.set_size(n,d);
    param_all.fill(0);
    return 0;
  } else {
    std::ifstream infile(f_init_params, std::ios_base::in);
    param_all.read_ascii(infile);
    assert(param_all.cols() == d);
    assert(param_all.rows() == n);
    return 1;
  }
}

void gmmreg_cpd::perform_transform(const vnl_vector<double> &x) {
  transformed_model = model + basis*param_all;
  //transformed_model = basis*param_all;
}

double gmmreg_cpd::bending_energy() {
  return vnl_trace(param_all.transpose() * kernel * param_all);
}

void gmmreg_cpd::save_results(const char* f_config,
    const vnl_vector<double>& params) {
  char f_transformed[256] = {0};
  GetPrivateProfileString(common_section, "transformed_model", NULL,
      f_transformed, 255, f_config);
  save_transformed(f_transformed, params, f_config );

  char f_final_params[256] = {0};
  GetPrivateProfileString(common_section, "final_params", NULL,
      f_final_params, 255, f_config);
  save_matrix(f_final_params, param_all);
}

void gmmreg_cpd::prepare_own_options(const char* f_config) {
  char s_EMtol[60] = {0}, s_anneal[60] = {0}, s_beta[60] = {0};
  char s_lambda[60] = {0}, s_outliers[60] = {0}, s_sigma[60] = {0};
  char s_tol[60] = {0}, s_viz[60] = {0};
  GetPrivateProfileString(section, "emtol", "1e-3", s_EMtol, 60, f_config);
  EMtol = atof(s_EMtol);
  GetPrivateProfileString(section, "anneal", "0.97", s_anneal, 60, f_config);
  anneal = atof(s_anneal);
  GetPrivateProfileString(section, "beta", "1", s_beta, 60, f_config);
  beta = atof(s_beta);
  GetPrivateProfileString(section, "lambda", "1", s_lambda, 60, f_config);
  lambda = atof(s_lambda);
  GetPrivateProfileString(section, "outliers", "0", s_outliers, 60, f_config);
  outliers = atoi(s_outliers);
  GetPrivateProfileString(section, "sigma", "1", s_sigma, 60, f_config);
  sigma = atof(s_sigma);
  GetPrivateProfileString(section, "tol", "1e-5", s_tol, 60, f_config);
  tol = atof(s_tol);
  max_iter = GetPrivateProfileInt(section, "max_iter", 150, f_config);
  max_em_iter = GetPrivateProfileInt(section, "max_em_iter", 150, f_config);
}
