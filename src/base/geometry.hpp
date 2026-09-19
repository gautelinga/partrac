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

inline double norm(const double x, const double y, const double z){
  return sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
}

inline double norm(const Vector3d &r){
  return r.norm();
}

inline double dist(const Uint i1, const Uint i2, std::vector<Vector3d>& x_rw){
  Vector3d dr = x_rw[i1]-x_rw[i2];
  return dr.norm();
}

inline double dist(const Vector3d &pta, const Vector3d &ptb){
  Vector3d dr = pta-ptb;
  return dr.norm();
}

inline double dot(const Vector3d &a, const Vector3d &b){
  return a.dot(b);
}

inline Vector3d diff(const Vector3d &a, const Vector3d &b){
  return a-b;
}

inline Vector3d cross(const Vector3d &a, const Vector3d &b){
  return a.cross(b);
}

inline double get_abs_angle(const Vector3d &a, const Vector3d &b){
  double costheta = a.dot(b)/(a.norm()*b.norm());
  return acos(costheta);
}

inline long double area(const Uint iedge, const Uint jedge,
                 std::vector<Vector3d>& x_rw,
                 const EdgesType& edges){
  Vector3d a = x_rw[edges[iedge].first[0]]-x_rw[edges[iedge].first[1]];
  Vector3d b = x_rw[edges[jedge].first[0]]-x_rw[edges[jedge].first[1]];
  return a.cross(b).norm()/2;
}

inline long double area(const Uint iface, std::vector<Vector3d>& x_rw,
                 const FacesType& faces, const EdgesType& edges){
  // To be decommissioned?
  Uint iedge = faces[iface].first[0];
  Uint jedge = faces[iface].first[1];
  return area(iedge, jedge, x_rw, edges);
}

inline std::vector<size_t> argsort_descending(const std::vector<double> &v){
  std::vector<size_t> idx(v.size());
  iota(idx.begin(), idx.end(), 0);
  stable_sort(idx.begin(), idx.end(),
              [&v](size_t i1, size_t i2) {return v[i1] > v[i2]; });
  return idx;
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
  //exit(1);
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

inline Vector3d vec_repl(const Uint inode,
                  std::vector<Vector3d>& x_rw,
                  const std::set<Uint> repl_nodes,
                  const Vector3d &x){
  if (contains(repl_nodes, inode))
    return x;
  return x_rw[inode];
}

inline Vector3d get_normal(const Uint iface,
                    const FacesType &faces,
                    const EdgesType &edges,
                    std::vector<Vector3d>& x_rw){
  Uint iedge = faces[iface].first[0];
  Uint jedge = faces[iface].first[1];
  Uint i00 = edges[iedge].first[0];
  Uint i01 = edges[iedge].first[1];
  Uint i10 = edges[jedge].first[0];
  Uint i11 = edges[jedge].first[1];

  Vector3d a = x_rw[i01]-x_rw[i00];
  Vector3d b = x_rw[i11]-x_rw[i10];
  Vector3d crossprod = a.cross(b);
  return crossprod/crossprod.norm();
}

inline Vector3d get_normal(const Uint jedge, const Uint kedge,
                    const EdgesType &edges,
                    std::vector<Vector3d>& x_rw,
                    const std::set<Uint> repl_nodes,
                    const Vector3d &x){
  Vector3d drj = vec_repl(edges[jedge].first[1], x_rw, repl_nodes, x)
    - vec_repl(edges[jedge].first[0], x_rw, repl_nodes, x);
  Vector3d drk = vec_repl(edges[kedge].first[1], x_rw, repl_nodes, x)
    - vec_repl(edges[kedge].first[0], x_rw, repl_nodes, x);
  return drj.cross(drk);
}

#endif
