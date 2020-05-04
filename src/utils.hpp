#include <iostream>
#include <math.h>

#ifndef __UTILS_HPP
#define __UTILS_HPP

using namespace std;

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

int imodulo(const int a, const int b) {
  return ((a % b) + b) % b;
}


double interpolate(const double x,
		   const double y,
		   const double z,
		   double*** C,
		   const int nx,
		   const int ny,
		   const int nz,
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
		    const int ind[3][2],
		    const double w[2][2][2]){
  double f = 0.0;
  for (int q0=0; q0<2; ++q0){
    for (int q1=0; q1<2; ++q1){
      for (int q2=0; q2<2; ++q2){
	f += C[ind[0][q0]][ind[1][q1]][ind[2][q2]]*w[q0][q1][q2];
      }
    }
  }
  return f;
}

#endif
