#ifndef GARS_UTIL
#define GARS_UTIL

#include <cmath>
#include <type_traits>


template <typename IntTy>
  requires std::is_integral_v<IntTy>
inline IntTy int_sqrt(IntTy x) {
  auto res = std::sqrt(static_cast<long double>(x));
  return static_cast<IntTy>(std::round(res));
}

#endif