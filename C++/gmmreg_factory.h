#ifndef GMMREG_FACTORY_H_
#define GMMREG_FACTORY_H_

#include <memory>
#include <string>

#include "gmmreg_base.h"

namespace gmmreg {

class GmmregFactory {
 public:
    static std::unique_ptr<Base> CreateInstance(const std::string& name);
};

}  // namespace gmmreg

#endif  // GMMREG_FACTORY_H_
