#ifndef gmmreg_rigid_h
#define gmmreg_rigid_h

#include <assert.h>

#include <iostream>
#include <fstream>

#include "gmmreg_base.h"
#include "gmmreg_rigid_func.h"

class gmmreg_rigid: public gmmreg_base {
 public:
  gmmreg_rigid() {
    strcpy(section, "GMMREG_OPT");
    func = new gmmreg_rigid_func;
  }
  ~gmmreg_rigid(){
    delete func;
  }

 private:
  vnl_vector<double> param_rigid;
  gmmreg_rigid_func *func;
  vnl_vector<long> nbd;
  vnl_vector<double> lower_bound, upper_bound;

  void set_bound();
  void start_registration(vnl_vector<double>&);
  int set_init_rigid(const char*);
  void set_param(vnl_vector<double>&);
  int set_init_params(const char*);
  void save_results(const char*, const vnl_vector<double>&);
  void prepare_basis_kernel() {};
  void prepare_param_gradient(bool) {};
  void perform_transform(const vnl_vector<double>&);
  double bending_energy() {
    return 0;
  };
  void compute_gradient(double lambda, const vnl_matrix<double>& gradient,
      vnl_matrix<double>& grad_all) {};
  void prepare_own_options(const char*);
};

#endif  // #ifndef gmmreg_rigid_h
