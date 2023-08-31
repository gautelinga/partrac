#ifndef __UTILS_HPP
#define __UTILS_HPP

#include <iostream>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <vector>
#include <list>
#include <set>
#include <Eigen/Dense>
#include <random>
#include <fstream>
#include "typedefs.hpp"
//#include "Interpol.hpp"
//#include "Parameters.hpp"

class Stamp {
public:
  Stamp(const double t_in, const std::string filename_in) : t(t_in), filename(filename_in) {};
  ~Stamp() {};
  double t;
  std::string filename;
};

class StampPair {
public:
  StampPair(const double t_prev, const std::string filename_prev,
	    const double t_next, const std::string filename_next) :
    prev(t_prev, filename_prev), next(t_next, filename_next) {};
  Stamp prev;
  Stamp next;
  double weight_next(const double t){ return (t-this->prev.t)/(this->next.t-this->prev.t); };
  double weight_prev(const double t){ return 1.-this->weight_next(t); };
};

static double modulox(const double x, const double L){
  if (x > 0){
    return fmod(x, L);
  }
  else {
    return fmod(x, L)+L;
  }
}

static Uint imodulo(const int a, const int b) {
  return ((a % b) + b) % b;
}

static double norm(const double x, const double y, const double z){
  return sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
}

static double norm(const Vector3d &r){
  return r.norm();
}

static double dist(const Uint i1, const Uint i2, std::vector<Vector3d>& x_rw){
  Vector3d dr = x_rw[i1]-x_rw[i2];
  return dr.norm();
}

static double dist(const Vector3d &pta, const Vector3d &ptb){
  Vector3d dr = pta-ptb;
  return dr.norm();
}

static double dot(const Vector3d &a, const Vector3d &b){
  return a.dot(b);
}

static Vector3d diff(const Vector3d &a, const Vector3d &b){
  return a-b;
}

static Vector3d cross(const Vector3d &a, const Vector3d &b){
  return a.cross(b);
}

static double get_abs_angle(const Vector3d &a, const Vector3d &b){
  double costheta = a.dot(b)/(a.norm()*b.norm());
  return acos(costheta);
}

static long double area(const Uint iedge, const Uint jedge,
                 std::vector<Vector3d>& x_rw,
                 const EdgesType& edges){
  Vector3d a = x_rw[edges[iedge].first[0]]-x_rw[edges[iedge].first[1]];
  Vector3d b = x_rw[edges[jedge].first[0]]-x_rw[edges[jedge].first[1]];
  return a.cross(b).norm()/2;
}

static long double area(const Uint iface, std::vector<Vector3d>& x_rw,
                 const FacesType& faces, const EdgesType& edges){
  // To be decommissioned?
  Uint iedge = faces[iface].first[0];
  Uint jedge = faces[iface].first[1];
  return area(iedge, jedge, x_rw, edges);
}

static std::vector<size_t> argsort_descending(const std::vector<double> &v){
  std::vector<size_t> idx(v.size());
  iota(idx.begin(), idx.end(), 0);
  stable_sort(idx.begin(), idx.end(),
              [&v](size_t i1, size_t i2) {return v[i1] > v[i2]; });
  return idx;
}

static Uint get_intersection(const std::array<Uint, 2> &a, const std::array<Uint, 2> &b){
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
  exit(0);
  return 0;
}

static Uint get_other(const Uint i, const Uint j, const Uint k){
  if (i==k)
    return j;
  assert(j==k);
  return i;
}

static double circumcenter(const Vector3d &A, const Vector3d &B, const Vector3d &C){
  Vector3d D((B-A).cross(C-A));
  double b = (A-C).norm();
  double c = (A-B).norm();
  double a = (B-C).norm();
  return 0.5*a*b*c/D.norm();
}

template<typename T>
bool contains(const std::set<T>& container, const T &elem){
  return container.find(elem) != container.end();
}

template<typename T1, typename T2>
bool contains(const std::map<T1, T2>& container, const T1 &elem){
  return container.find(elem) != container.end();
}

static Vector3d vec_repl(const Uint inode,
                  std::vector<Vector3d>& x_rw,
                  const std::set<Uint> repl_nodes,
                  const Vector3d &x){
  if (contains(repl_nodes, inode))
    return x;
  return x_rw[inode];
}

static Vector3d get_normal(const Uint iface,
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

static Vector3d get_normal(const Uint jedge, const Uint kedge,
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

static double getd(std::map<std::string, std::string> &expr_params, const std::string key){
  if (expr_params.find(key) != expr_params.end()){
    return stod(expr_params[key]);
  }
  else {
    std::cout << "Missing key: " << key << std::endl;
    exit(0);
    return 0.;
  }
}

static double geti(std::map<std::string, std::string> &expr_params, const std::string key){
  if (expr_params.find(key) != expr_params.end()){
    return stoi(expr_params[key]);
  }
  else {
    std::cout << "Missing key: " << key << std::endl;
    exit(0);
    return 0;
  }
}

static std::vector<std::string> split_string(const std::string s, const std::string delim){
  std::vector<std::string> s_;
  auto start = 0U;
  auto end = s.find(delim);
  while (end != std::string::npos){
    s_.push_back(s.substr(start, end-start));
    start = end + delim.length();
    end = s.find(delim, start);
  }
  s_.push_back(s.substr(start, end));
  return s_;
}

static bool contains(const std::string s, const std::string c){
  return (s.find(c) != std::string::npos);
}

static std::vector<double> getdvec(std::map<std::string, std::string> &expr_params,
                       const std::string key){
  std::string token = ",";
  std::vector<double> dvec;
  if (expr_params.find(key) != expr_params.end()){
    std::string str = expr_params[key];
    while (str.size()){
      Uint index = str.find(token);
      if (index != std::string::npos){
        dvec.push_back(stod(str.substr(0, index)));
        str = str.substr(index+token.size());
        if (str.size() == 0) {
          dvec.push_back(stod(str));
        }
      }
      else {
        dvec.push_back(stod(str));
        str = "";
      }
    }
  }
  else {
    std::cout << "Missing key: " << key << std::endl;
    exit(0);
  }
  return dvec;
}

static std::vector<int> getivec(std::map<std::string, std::string> &expr_params,
                    const std::string key){
  std::string token = ",";
  std::vector<int> ivec;
  if (expr_params.find(key) != expr_params.end()){
    std::string str = expr_params[key];
    while (str.size()){
      Uint index = str.find(token);
      if (index != std::string::npos){
        ivec.push_back(stoi(str.substr(0, index)));
        str = str.substr(index+token.size());
        if (str.size() == 0) {
          ivec.push_back(stoi(str));
        }
      }
      else {
        ivec.push_back(stoi(str));
        str = "";
      }
    }
  }
  else {
    std::cout << "Missing key: " << key << std::endl;
    exit(0);
  }
  return ivec;
}

class PointValues {
public:
  PointValues(const double U0) : U0(U0) {
    U = {0., 0., 0.};
    A = {0., 0., 0.};
    gradU << 0., 0., 0., 0., 0., 0., 0., 0., 0.; 
    gradA << 0., 0., 0., 0., 0., 0., 0., 0., 0.; 
    P = 0.;
  };
  Vector3d U;
  Vector3d A;
  Matrix3d gradU;
  Matrix3d gradA;
  double P;
  Vector3d get_u() { return U0 * U; };
  Vector3d get_Ju() { return U0 * U0 * gradU * U; }; // check
  Vector3d get_a() { return U0 * A; };
private:
  double U0;
};

#endif