#ifndef __INTERVALS_HPP
#define __INTERVALS_HPP

#include <limits>
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

#endif
