#ifndef __GEOMETRY_HPP
#define __GEOMETRY_HPP

#include <cmath>
#include <numeric>
#include <algorithm>
#include <array>
#include <cassert>
#include <iostream>
#include <set>
#include <vector>
#include "typedefs.hpp"
#include "strings.hpp"

inline double modulox(const double x, const double L){
  if (x > 0){
    return fmod(x, L);
  }
  else {
    return fmod(x, L)+L;
  }
}

inline Uint imodulo(const int a, const int b) {
  // In range: no division
  if (static_cast<unsigned>(a) < static_cast<unsigned>(b)) return a;
  return ((a % b) + b) % b;
}

inline double dist(const Vector3d &pta, const Vector3d &ptb){
  Vector3d dr = pta-ptb;
  return dr.norm();
}

inline Uint get_intersection(const std::array<Uint, 2> &a, const std::array<Uint, 2> &b){
  for (std::array<Uint, 2>::const_iterator ait=a.begin();
       ait != a.end(); ++ait){
    for (std::array<Uint, 2>::const_iterator bit=b.begin();
         bit != b.end(); ++bit){
      if (*ait == *bit){
        return *ait;
      }
    }
  }
  std::cout << "Error: found no intersection." << std::endl;
  // GL: Hack to avoid crashing. This function is only used for curvature calculations!
  return 0;
}

inline Uint get_other(const Uint i, const Uint j, const Uint k){
  if (i==k)
    return j;
  assert(j==k);
  return i;
}

inline double circumcenter(const Vector3d &A, const Vector3d &B, const Vector3d &C){
  Vector3d D((B-A).cross(C-A));
  double b = (A-C).norm();
  double c = (A-B).norm();
  double a = (B-C).norm();
  return 0.5*a*b*c/D.norm();
}

#endif
