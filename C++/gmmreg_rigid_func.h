#ifndef gmmreg_rigid_func_h
#define gmmreg_rigid_func_h

#include <vnl/vnl_cost_function.h>
#include "gmmreg_base.h"

class gmmreg_rigid_func : public vnl_cost_function {
 public:
  gmmreg_rigid_func(): vnl_cost_function() {}

  double eval(double &, vnl_matrix<double> &);
  double f(const vnl_vector<double>& x);
  void gradf(const vnl_vector<double>& x, vnl_vector<double>& g);

  inline void set_scale(double scale) {
    this->scale = scale;
  }
  inline double get_scale() const {
    return this->scale;
  }

  gmmreg_base* gmmreg;
  inline void set_gmmreg(gmmreg_base* gmmreg) {
    this->gmmreg = gmmreg;
    this->m = gmmreg->m;
    this->d = gmmreg->d;
    if (d == 2) {
      set_number_of_unknowns(3);
    } else if (d == 3) {
      set_number_of_unknowns(7);
    }
    gradient.set_size(m, d);
  }

  ~gmmreg_rigid_func() {}

 protected:
  vnl_matrix<double> gradient;

 private:
  double scale, lambda;
  int m, d;
};

#endif  // #ifndef gmmreg_rigid_func_h_
