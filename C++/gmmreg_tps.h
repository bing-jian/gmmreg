#ifndef GMMREG_TPS_H_
#define GMMREG_TPS_H_

#include <vector>

#include "gmmreg_base.h"
#include "gmmreg_tps_func.h"

namespace gmmreg {

class TpsRegistration: public Base {
 public:
  TpsRegistration() {
    strcpy(section_, "GMMREG_OPT");
  }
  virtual ~TpsRegistration() {
    delete func_;
  }

  // Set initial TPS params directly (no file needed).
  // affine: (d+1) x d;  tps: (n-d-1) x d.
  // Call after Prepare() so that n_ and d_ are known.
  void SetInitParams(const vnl_matrix<double>& affine,
                     const vnl_matrix<double>& tps) {
    param_affine_ = affine;
    param_tps_    = tps;
    param_all_.set_size(n_, d_);
    param_all_.update(param_affine_);
    param_all_.update(param_tps_, d_ + 1);
  }

 protected:
  ThinPlateSplineFunc* func_;

 private:
  vnl_matrix<double> param_affine_, param_tps_, param_all_;
#ifdef PRE_COMPUTE_TPS_BASIS_AND_KERNEL
  vnl_matrix<double> basis_;
#endif
  // Kept intermediate matrices in case we need to optimize
  // matrix chain multiplication.
  vnl_matrix<double> K_, U_, PP_, PP_transpose_, U_transpose_;

  std::vector<float> v_lambda_;
  std::vector<int> v_affine_;

  void StartRegistration(vnl_vector<double>&) override;
  void ApplyInitParams(const RegistrationInput&) override;
  int SetInitAffine(const char* filename);
  int SetInitTps(const char* filename);
  void SetParam(vnl_vector<double>& x0);
  void SetAffineAndTps(const vnl_vector<double>&);
  int SetInitParams(const char* filename) override;
  void SaveResults(const char* f_config, const vnl_vector<double>&) override;

  int PrepareInput(const char* input_config) override;
  void PrepareBasisKernel() override;
  void PrepareParamGradient(bool);
  void PerformTransform(const vnl_vector<double>&) override;
  double BendingEnergy() override;
  void ComputeGradient(const double lambda, const vnl_matrix<double>& gradient,
      vnl_matrix<double>& grad_all) override;
  void PrepareOwnOptions(const char* f_config) override;
};

class TpsRegistration_L2: public TpsRegistration {
 public:
  TpsRegistration_L2(): TpsRegistration() {
    func_ = new ThinPlateSplineFunc_L2;
  }
};

class TpsRegistration_KC: public TpsRegistration {
 public:
  TpsRegistration_KC(): TpsRegistration() {
    func_ = new ThinPlateSplineFunc_KC;
  }
};

}  // namespace gmmreg
#endif  // GMMREG_TPS_H_
