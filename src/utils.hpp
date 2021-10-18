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
#include "Interpol.hpp"

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

static double interpolate(const double x,
		   const double y,
		   const double z,
		   double*** C,
		   const Uint nx,
		   const Uint ny,
		   const Uint nz,
		   const double Lx,
		   const double Ly,
		   const double Lz,
		   const bool verbose
		   ){
  int ix_lo, ix_hi, iy_lo, iy_hi, iz_lo, iz_hi;
  double wx, wy, wz;
  double f000, f001, f010, f011, f100, f101, f110, f111;
  double f;
  double dx = Lx/nx;
  double dy = Ly/ny;
  double dz = Lz/nz;

  int ix_est = floor(x/dx);
  int iy_est = floor(y/dy);
  int iz_est = floor(z/dz);

  ix_lo = imodulo(ix_est, nx);
  ix_hi = imodulo(ix_lo + 1, nx);
  iy_lo = imodulo(iy_est, ny);
  iy_hi = imodulo(iy_lo + 1, ny);
  iz_lo = imodulo(iz_est, nz);
  iz_hi = imodulo(iz_lo + 1, nz);

  wx = (x-dx*ix_est)/dx;
  wy = (y-dy*iy_est)/dy;
  wz = (z-dz*iz_est)/dz;

  if (verbose)
    std::cout << wx << " " << wy << " " << wz << std::endl;

  // Nodal values
  f000 = C[ix_lo][iy_lo][iz_lo];
  f001 = C[ix_lo][iy_lo][iz_hi];
  f010 = C[ix_lo][iy_hi][iz_lo];
  f011 = C[ix_lo][iy_hi][iz_hi];
  f100 = C[ix_hi][iy_lo][iz_lo];
  f101 = C[ix_hi][iy_lo][iz_hi];
  f110 = C[ix_hi][iy_hi][iz_lo];
  f111 = C[ix_hi][iy_hi][iz_hi];

  // Trilinear interpolaton
  f = f000*(1-wx)*(1-wy)*(1-wz)
    + f001*(1-wx)*(1-wy)*  wz
    + f010*(1-wx)*  wy  *(1-wz)
    + f011*(1-wx)*  wy  *  wz
    + f100*  wx  *(1-wy)*(1-wz)
    + f101*  wx  *(1-wy)*  wz
    + f110*  wx  *  wy  *(1-wz)
    + f111*  wx  *  wy  *  wz;

  return f;
}

static double weighted_sum(double*** C,
                    const Uint ind[3][2],
                    const double w[2][2][2]){
  double f = 0.0;
  for (Uint q0=0; q0<2; ++q0){
    for (Uint q1=0; q1<2; ++q1){
      for (Uint q2=0; q2<2; ++q2){
        f += C[ind[0][q0]][ind[1][q1]][ind[2][q2]]*w[q0][q1][q2];
      }
    }
  }
  return f;
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

static void test_interpolation(Uint num_points, Interpol *intp,
                        const std::string &newfolder, const double t0,
                        std::mt19937 &gen){
  Uint n = 0;
  Uint n_inside = 0;

  Vector3d x_min = intp->get_x_min();
  Vector3d x_max = intp->get_x_max();

  double Lx = intp->get_Lx();
  double Ly = intp->get_Ly();
  double Lz = intp->get_Lz();

  std::cout << "Lx=" << Lx << ", Ly=" << Ly << ", Lz=" << Lz << std::endl;

  intp->update(t0);

  // std::ofstream nodalfile(newfolder + "/nodal_values.dat");
  // for (Uint ix=0; ix<intp->get_nx(); ++ix){
  //   for (Uint iy=0; iy<intp->get_ny(); ++iy){
  //     for (Uint iz=0; iz<intp->get_nz(); ++iz){
  //       bool inside = intp->get_nodal_inside(ix, iy, iz);
  //       Vector3d u(intp->get_nodal_ux(ix, iy, iz),
  //                  intp->get_nodal_uy(ix, iy, iz),
  //                  intp->get_nodal_uz(ix, iy, iz));
  //       nodalfile << ix << " " << iy << " " << iz << " " << inside << " "
  //                 << u[0] << " " << u[1] << " " << u[2] << std::endl;
  //     }
  //   }
  // }
  // nodalfile.close();

  std::uniform_real_distribution<> uni_dist_x(x_min[0], x_max[0]);
  std::uniform_real_distribution<> uni_dist_y(x_min[1], x_max[1]);
  std::uniform_real_distribution<> uni_dist_z(x_min[2], x_max[2]);

  std::ofstream ofile(newfolder + "/interpolation.txt");
  std::string sep = ",";
  ofile << "x" << sep << "y" << sep << "z" << sep
        << "ux" << sep << "uy" << sep << "uz" << sep
        << "rho" << sep << "p" << sep << "divu" << sep
        << "vortz" << sep
        << "uxx" << sep << "uxy" << sep << "uxz" << sep
        << "uyx" << sep << "uyy" << sep << "uyz" << sep
        << "uzx" << sep << "uzy" << sep << "uzz"
        << std::endl;
  while (n < num_points){
    Vector3d x(uni_dist_x(gen), uni_dist_y(gen), uni_dist_z(gen));
    //std::cout << x << std::endl;
    intp->probe(x);
    if (intp->inside_domain()){
      Vector3d u = intp->get_u();
      double rho = intp->get_rho();
      double p = intp->get_p();
      double divu = intp->get_divu();
      double vortz = intp->get_vortz();
      ofile << x[0] << sep << x[1] << sep << x[2] << sep
            << u[0] << sep << u[1] << sep << u[2] << sep
            << rho << sep << p << sep << divu << sep
            << vortz << sep
            << intp->get_uxx() << sep << intp->get_uxy() << sep << intp->get_uxz() << sep
            << intp->get_uyx() << sep << intp->get_uyy() << sep << intp->get_uyz() << sep
            << intp->get_uzx() << sep << intp->get_uzy() << sep << intp->get_uzz()
            << std::endl;
      ++n_inside;
    }
    ++n;
  }
  std::cout << "Inside: " << n_inside << "/" << n << std::endl;
  std::cout << "Approximate volume: " << (n_inside*Lx*Ly*Lz)/n << std::endl;
  std::cout << "Approximate area:   " << (n_inside*Lx*Ly)/n << std::endl;
  ofile.close();
}

#endif
