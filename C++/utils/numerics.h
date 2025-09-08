#ifndef UTILS_NUMERICS_H_
#define UTILS_NUMERICS_H_

namespace gmmreg {

template<typename T>
inline T exp_256(T x) {
  x = 1.0 + x / 256.0;
  x *= x; x *= x; x *= x; x *= x;
  x *= x; x *= x; x *= x; x *= x;
  return x;
}

}  // namespace gmmreg

#endif  // UTILS_NUMERICS_H_
