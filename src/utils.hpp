#include <iostream>
#include <math.h>
#include <numeric>
#include <algorithm>
#include <vector>
#include <list>
#include <set>
#include <Eigen/Dense>

#ifndef __UTILS_HPP
#define __UTILS_HPP

using namespace std;

typedef size_t Uint;

typedef vector<pair<array<Uint, 2>, double>> EdgesType;
typedef vector<pair<array<Uint, 3>, double>> FacesType;
typedef list<Uint> FacesListType;
typedef list<Uint> EdgesListType;
typedef vector<FacesListType> Edge2FacesType;
typedef vector<EdgesListType> Node2EdgesType;
typedef vector<map<Uint, double>> InteriorAnglesType;
typedef Eigen::Vector3d Vector3d;
typedef Eigen::Matrix3d Matrix3d;

class Stamp {
public:
  Stamp(const double t_in, const string filename_in) : t(t_in), filename(filename_in) {};
  ~Stamp() {};
  double t;
  string filename;
};

class StampPair {
public:
  StampPair(const double t_prev, const string filename_prev,
	    const double t_next, const string filename_next) :
    prev(t_prev, filename_prev), next(t_next, filename_next) {};
  Stamp prev;
  Stamp next;
  double weight_next(const double t){ return (t-this->prev.t)/(this->next.t-this->prev.t); };
  double weight_prev(const double t){ return 1.-this->weight_next(t); };
};

double modulox(const double x, const double L){
  if (x > 0){
    return fmod(x, L);
  }
  else {
    return fmod(x, L)+L;
  }
}

Uint imodulo(const int a, const int b) {
  return ((a % b) + b) % b;
}

double interpolate(const double x,
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

double weighted_sum(double*** C,
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

double norm(const double x, const double y, const double z){
  return sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
}

double norm(const Vector3d &r){
  return r.norm();
}

double dist(const Uint i1, const Uint i2, Vector3d* x_rw){
  Vector3d dr = x_rw[i1]-x_rw[i2];
  return dr.norm();
}

double dist(const Vector3d &pta, const Vector3d &ptb){
  Vector3d dr = pta-ptb;
  return dr.norm();
}

double dot(const Vector3d &a, const Vector3d &b){
  return a.dot(b);
}

Vector3d diff(const Vector3d &a, const Vector3d &b){
  return a-b;
}

Vector3d cross(const Vector3d &a, const Vector3d &b){
  return a.cross(b);
}

double get_abs_angle(const Vector3d &a, const Vector3d &b){
  double costheta = a.dot(b)/(a.norm()*b.norm());
  return acos(costheta);
}

long double area(const Uint iedge, const Uint jedge,
                 Vector3d* x_rw,
                 const EdgesType& edges){
  Vector3d a = x_rw[edges[iedge].first[0]]-x_rw[edges[iedge].first[1]];
  Vector3d b = x_rw[edges[jedge].first[0]]-x_rw[edges[jedge].first[1]];
  return a.cross(b).norm()/2;
}

long double area(const Uint iface, Vector3d* x_rw,
                 const FacesType& faces, const EdgesType& edges){
  Uint iedge = faces[iface].first[0];
  Uint jedge = faces[iface].first[1];
  return area(iedge, jedge, x_rw, edges);
}

vector<size_t> argsort_descending(const vector<double> &v){
  vector<size_t> idx(v.size());
  iota(idx.begin(), idx.end(), 0);
  stable_sort(idx.begin(), idx.end(),
              [&v](size_t i1, size_t i2) {return v[i1] > v[i2]; });
  return idx;
}

Uint get_intersection(const array<Uint, 2> &a, const array<Uint, 2> &b){
  for (array<Uint, 2>::const_iterator ait=a.begin();
       ait != a.end(); ++ait){
    for (array<Uint, 2>::const_iterator bit=b.begin();
         bit != b.end(); ++bit){
      if (*ait == *bit){
        return *ait;
      }
    }
  }
  cout << "Error: found no intersection." << endl;
  exit(0);
  return 0;
}

Uint get_other(const Uint i, const Uint j, const Uint k){
  if (i==k)
    return j;
  assert(j==k);
  return i;
}

double circumcenter(const Vector3d &A, const Vector3d &B, const Vector3d &C){
  Vector3d D((B-A).cross(C-A));
  double b = (A-C).norm();
  double c = (A-B).norm();
  double a = (B-C).norm();
  return 0.5*a*b*c/D.norm();
}

template<typename T>
bool contains(const set<T>& container, const T &elem){
  return container.find(elem) != container.end();
}

template<typename T1, typename T2>
bool contains(const map<T1, T2>& container, const T1 &elem){
  return container.find(elem) != container.end();
}

Vector3d vec_repl(const Uint inode,
                  Vector3d* x_rw,
                  const set<Uint> repl_nodes,
                  const Vector3d &x){
  if (contains(repl_nodes, inode))
    return x;
  return x_rw[inode];
}

Vector3d get_normal(const Uint iface,
                    const FacesType &faces,
                    const EdgesType &edges,
                    Vector3d* x_rw){
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

Vector3d get_normal(const Uint jedge, const Uint kedge,
                    const EdgesType &edges,
                    Vector3d* x_rw,
                    const set<Uint> repl_nodes,
                    const Vector3d &x){
  Vector3d drj = vec_repl(edges[jedge].first[1], x_rw, repl_nodes, x)
    - vec_repl(edges[jedge].first[0], x_rw, repl_nodes, x);
  Vector3d drk = vec_repl(edges[kedge].first[1], x_rw, repl_nodes, x)
    - vec_repl(edges[kedge].first[0], x_rw, repl_nodes, x);
  return drj.cross(drk);
}

#endif
