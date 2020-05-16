#include <iostream>
#include <math.h>
#include <numeric>
#include <algorithm>
#include <vector>
#include <list>

#ifndef __UTILS_HPP
#define __UTILS_HPP

using namespace std;

typedef size_t Uint;

typedef vector<pair<array<Uint, 2>, double>> EdgesType;
typedef vector<pair<array<Uint, 3>, double>> FacesType;
typedef list<Uint> FacesListType;
typedef vector<FacesListType> Edge2FacesType;

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

double norm(const array<double, 3> r){
  return norm(r[0], r[1], r[2]);
}

double dist(const Uint i1, const Uint i2, double* x_rw, double* y_rw, double* z_rw){
  return norm(x_rw[i1]-x_rw[i2], y_rw[i1]-y_rw[i2], z_rw[i1]-z_rw[i2]);
}

double dist(const array<double, 3> pta, const array<double, 3> ptb){
  double x_a = pta[0];
  double y_a = pta[1];
  double z_a = pta[2];
  double x_b = ptb[0];
  double y_b = ptb[1];
  double z_b = ptb[2];
  return norm(x_a-x_b, y_a-y_b, z_a-z_b);
}

array<double, 3> cross(const array<double, 3> a, const array<double, 3> b){
  return {a[1]*b[2]-a[2]*b[1],
          a[2]*b[0]-a[0]*b[2],
          a[0]*b[1]-a[1]*b[0]};
}

double area(const Uint iedge, const Uint jedge, double* x_rw, double* y_rw, double* z_rw, const EdgesType& edges){
  array<double, 3> a = {x_rw[edges[iedge].first[0]]-x_rw[edges[iedge].first[1]],
                        y_rw[edges[iedge].first[0]]-y_rw[edges[iedge].first[1]],
                        z_rw[edges[iedge].first[0]]-z_rw[edges[iedge].first[1]]};
  array<double, 3> b = {x_rw[edges[jedge].first[0]]-x_rw[edges[jedge].first[1]],
                        y_rw[edges[jedge].first[0]]-y_rw[edges[jedge].first[1]],
                        z_rw[edges[jedge].first[0]]-z_rw[edges[jedge].first[1]]};
  return norm(cross(a, b))/2;
}

double area(const Uint iface, double* x_rw, double* y_rw, double* z_rw,
            const FacesType& faces, const EdgesType& edges){
  Uint iedge = faces[iface].first[0];
  Uint jedge = faces[iface].first[1];
  return area(iedge, jedge, x_rw, y_rw, z_rw, edges);
}

vector<size_t> argsort_descending(const vector<double> &v){
  vector<size_t> idx(v.size());
  iota(idx.begin(), idx.end(), 0);
  stable_sort(idx.begin(), idx.end(),
              [&v](size_t i1, size_t i2) {return v[i1] > v[i2]; });
  return idx;
}

#endif
