#ifndef gmmreg_grbf_func_h
#define gmmreg_grbf_func_h

#include <vnl/vnl_cost_function.h>
#include "gmmreg_base.h"

class gmmreg_grbf_func : public vnl_cost_function {
 public:
  gmmreg_grbf_func(): vnl_cost_function(), m(0), n(0), d(0) {}

  virtual double eval(double& f1, double& f2,
      vnl_matrix<double>& g1, vnl_matrix<double>& g2) = 0;
  double f(const vnl_vector<double>& x);
  void gradf(const vnl_vector<double>& x, vnl_vector<double>& g);

  gmmreg_base* gmmreg;
  inline void set_gmmreg(gmmreg_base* gmmreg) {
    this->gmmreg = gmmreg;
    this->m = gmmreg->m;
    this->n = gmmreg->n;
    this->d = gmmreg->d;
    gradient1.set_size(m, d);
    gradient2.set_size(m, d);
    set_number_of_unknowns(n * d); //dim = n*d;
  }
  inline void set_scale(double scale) {
    this->scale = scale;
  }
  inline double get_scale() const {
    return this->scale;
  }
  inline void set_lambda(double lambda) {
    this->lambda = lambda;
  }
  inline double get_lambda() const {
    return this->lambda;
  }
  inline void set_beta(double beta) {
    this->beta = beta;
  }
  inline double get_beta() const {
    return this->beta;
  }
  void prepare_param_gradient();
  virtual ~gmmreg_grbf_func() {}

 protected:
  vnl_matrix<double> gradient;

 private:
  double scale, lambda, beta;
  int m, n, d;
  vnl_matrix<double> gradient1, gradient2, grad_all;
};

class gmmreg_grbf_L2_func : public gmmreg_grbf_func {
  double eval(double&f1, double &f2,
      vnl_matrix<double> &g1, vnl_matrix<double> &g2);
};

class gmmreg_grbf_KC_func : public gmmreg_grbf_func {
  double eval(double&f1, double &f2,
      vnl_matrix<double> &g1, vnl_matrix<double> &g2);
};

#endif //#ifndef gmmreg_grbf_func_h_
