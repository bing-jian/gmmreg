#ifndef gmmreg_base_h
#define gmmreg_base_h

#ifdef WIN32
#include <windows.h>
#else
#include "port_ini.h"
#endif

#include <vector>
#include <vnl/vnl_vector.h>
#include <vnl/vnl_matrix.h>

class gmmreg_base {
 public:
  gmmreg_base() {
    strcpy(common_section, "FILES");
  }
  int m, n, s, d;
  /* m: # of points in model */
  /* s: # of points in scene */
  /* n: # of points in ctrl_pts */
  /* d: dimensionality, e.g. 2 for 2D points, 3 for 3D points */

  // each row is a sample point
  vnl_matrix<double> model, scene, ctrl_pts, transformed_model;

  void run(const char* f_config);
  virtual void perform_transform(const vnl_vector<double>&) = 0;
  virtual double bending_energy() = 0; //serving as a regularization term
  virtual void compute_gradient(double lambda,
      const vnl_matrix<double>& gradient, vnl_matrix<double>& grad_all) = 0;

  virtual ~gmmreg_base() {}

 protected:
  double sigma, lambda;
  vnl_matrix<double> kernel;
  int b_normalize;
  vnl_vector<double> model_centroid, scene_centroid;
  char section[80], common_section[80];

  unsigned int level;
  std::vector<double> v_scale;
  std::vector<int> v_func_evals;

  // load input data from files
  virtual int prepare_input(const char* input_config);
  int set_ctrl_pts(const char* filename);
  void save_transformed(const char* filename,
      const vnl_vector<double>&, const char* f_config);
  void multi_scale_options(const char* f_config);
  void denormalize_all();
  int initialize(const char* f_config);

 private:
  double model_scale, scene_scale;
  void prepare_common_options(const char* f_config);
  virtual void prepare_own_options(const char* f_config) = 0;
  virtual void prepare_basis_kernel() = 0;
  virtual int set_init_params(const char* filename) = 0;
  virtual void save_results(const char* filename,
      const vnl_vector<double>&) = 0;
  virtual void start_registration(vnl_vector<double>& params) = 0;
};

#endif // #ifndef gmmreg_base_h
