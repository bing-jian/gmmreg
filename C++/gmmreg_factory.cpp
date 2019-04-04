#include "gmmreg_factory.h"

#include <memory>
#include <string>

#include "gmmreg_cpd.h"
#include "gmmreg_grbf.h"
#include "gmmreg_rigid.h"
#include "gmmreg_tps.h"

namespace gmmreg {

std::unique_ptr<Base> GmmregFactory::CreateInstance(const std::string& name) {
  Base * instance = nullptr;
  if (name == "em_tps") {
    instance = new CoherentPointDriftTps();
  } else if (name == "em_grbf") {
    instance = new CoherentPointDriftGrbf();
  } else if (name == "tps_l2") {
    instance = new TpsRegistration_L2();
  } else if (name == "tps_kc") {
    instance = new TpsRegistration_KC();
  } else if (name == "grbf_l2") {
    instance = new GrbfRegistration_L2();
  } else if (name == "grbf_kc") {
    instance = new GrbfRegistration_KC();
  } else if (name == "rigid") {
    instance = new RigidRegistration();
  }
  if (instance != nullptr) {
    return std::unique_ptr<Base>(instance);
  } else {
    return nullptr;
  }
}

}  // namespace gmmreg
