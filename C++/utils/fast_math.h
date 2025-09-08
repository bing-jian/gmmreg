#ifndef GMMREG_UTILS_FAST_MATH_H_
#define GMMREG_UTILS_FAST_MATH_H_

namespace gmmreg {
namespace utils {

float fastExp4(register float x);
double FastExp(double x);
double exp_256(double x);

}  // namespace utils
}  // namespace gmmreg

#endif  // GMMREG_UTILS_FAST_MATH_H_
