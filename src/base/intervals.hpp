#ifndef __INTERVALS_HPP
#define __INTERVALS_HPP

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <string>
#include "typedefs.hpp"

// Timesteps in an interval; at least 1, saturating
inline Uint steps_per(const double intv, const double dt){
  const double n = intv/dt;
  if (!(n > 1.)) return 1;  // shorter than a timestep, or NaN
  if (n >= double(std::numeric_limits<Uint>::max()))
    return std::numeric_limits<Uint>::max();
  return Uint(n);
}

// Whether step it is on the interval; off if intv <= 0
inline bool at_interval(const Uint it, const double intv, const double dt){
  if (!(intv > 0.)) return false;
  return it % steps_per(intv, dt) == 0;
}

// Time label: six decimals, more if dt needs them
inline std::string time_name(const double t, const double dt){
  const int digits = dt > 0. ? std::clamp(int(std::ceil(-std::log10(dt))) + 2, 6, 17) : 6;
  char s[64];
  std::snprintf(s, sizeof s, "%.*f", digits, t);
  return s;
}

#endif
