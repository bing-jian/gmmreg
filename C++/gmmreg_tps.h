#ifndef gmmreg_tps_h
#define gmmreg_tps_h

#include <vector>

#include "gmmreg_base.h"
#include "gmmreg_tps_func.h"

class gmmreg_tps: public gmmreg_base {
 public:
  gmmreg_tps() {
    strcpy(section, "GMMREG_OPT");
  }
  virtual ~gmmreg_tps() {
    delete func;
  }

 protected:
  gmmreg_tps_func *func;

 private:
  vnl_matrix<double> param_affine, param_tps;
  vnl_matrix<double> after_tps, basis, param_all;
  std::vector<double> v_lambda;
  std::vector<int> v_affine;

  void start_registration(vnl_vector<double>&);
  int set_init_affine(const char* filename);
  int set_init_tps(const char* filename);
  void set_param(vnl_vector<double>& x0);
  void set_affine_and_tps(const vnl_vector<double>&);
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

class gmmreg_tps_L2: public gmmreg_tps {
 public:
  gmmreg_tps_L2(): gmmreg_tps() {
    func = new gmmreg_tps_L2_func;
  }
};

class gmmreg_tps_KC: public gmmreg_tps {
 public:
  gmmreg_tps_KC(): gmmreg_tps() {
    func = new gmmreg_tps_KC_func;
  }
};

#endif //#ifndef gmmreg_tps_h
