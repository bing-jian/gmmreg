#ifndef gmmreg_grbf_h
#define gmmreg_grbf_h

#include <vector>

#include "gmmreg_base.h"
#include "gmmreg_grbf_func.h"

class gmmreg_grbf: public gmmreg_base {
 public:
  gmmreg_grbf() {
    strcpy(section, "GMMREG_OPT");
  }
  virtual ~gmmreg_grbf() {
    delete func;
  }

 protected:
  gmmreg_grbf_func *func;

 private:
  vnl_matrix<double> param_grbf;
  vnl_matrix<double> after_grbf, basis, param_all;
  std::vector<double> v_beta;
  double beta;
  std::vector<double> v_lambda;
  std::vector<int> v_affine;

  void start_registration(vnl_vector<double>&);
  int set_init_grbf(const char* filename);
  void set_param(vnl_vector<double>& x0);
  void set_grbf(const vnl_vector<double>&);
  int set_init_params(const char* filename);
  void save_results(const char* f_config, const vnl_vector<double>&);

  int prepare_input(const char* input_config);
  void prepare_basis_kernel();
  void prepare_param_gradient(bool);
  void perform_transform(const vnl_vector<double>&);
  double bending_energy();
  void compute_gradient(double lambda, const vnl_matrix<double>& gradient,
      vnl_matrix<double>& grad_all);
  void prepare_own_options(const char* f_config);
};

class gmmreg_grbf_L2: public gmmreg_grbf {
 public:
  gmmreg_grbf_L2(): gmmreg_grbf() {
    func = new gmmreg_grbf_L2_func;
  }
};

class gmmreg_grbf_KC: public gmmreg_grbf {
 public:
  gmmreg_grbf_KC(): gmmreg_grbf() {
    func = new gmmreg_grbf_KC_func;
  }
};

#endif  // #ifndef gmmreg_grbf_h
