#ifndef __CURV_LIMIT_HPP
#define __CURV_LIMIT_HPP

#include <algorithm>

// Edge lengths shrink where a strip bends: by 1/min(1 + factor kappa, cap);
// no limit with factor 0. Strips only: splits alone cannot grade a sheet
// that steeply (Topology::check_dim refuses)
struct CurvLimit {
  double factor = 0.;
  double cap = 16.;
  double scale(const double kappa) const { return 1./std::min(1. + factor*kappa, cap); }
};

#endif
