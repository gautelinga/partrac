#ifndef TET_HPP
#define TET_HPP

#include <array>
#include "typedefs.hpp"

#ifdef USE_DOLFIN
namespace dolfin { class Cell; }
#endif

class alignas(64) Tet
{

public:

  static constexpr int n_verts = 4;
  static constexpr std::size_t n_dofs_max = 10;   // quadbasis writes this many

  Tet() {}
  // The four vertex coordinate triples, in dolfin's local vertex order
  Tet(const double* x0, const double* x1, const double* x2, const double* x3);
#ifdef USE_DOLFIN
  Tet(const dolfin::Cell& cell);
#endif

  void xyz2bary(double x, double y, double z,
                double &r, double &s, double &t, double &u) const;

  // Writes the barycentrics, inside or not
  bool contains(const Vector3d& x, std::array<double, 4>& bary) const;

  // Gradient of barycentric k: inward normal of the facet opposite vertex k
  Vector3d bary_grad(const int k) const;

  void linearbasis(double r, double s, double t, double u,
                   double *N) const;

  void linearderiv(double, double, double, double,
                   double *Nx,
                   double *Ny,
                   double *Nz) const;

  void quadbasis(double r, double s, double t, double u,
                 double *N) const;

  void quadderiv(double r,double s,double t,double u,
                 double *Nx,
                 double *Ny,
                 double *Nz) const;

private:

  // Barycentric data first
  double x0_, y0_, z0_;
  double g2x_, g2y_, g2z_;
  double g3x_, g3y_, g3z_;
  double g4x_, g4y_, g4z_;
  double g1x_, g1y_, g1z_;

  static constexpr std::array<int, 10> perm_ = {-1, -1, -1, -1, 9, 6, 8, 7, 5, 4};

public:

  // quadbasis slots of the midpoints of edges 01, 02, 03, 12, 13, 23
  static constexpr std::array<int, 6> mid_ = {perm_[4], perm_[6], perm_[7], perm_[5], perm_[8], perm_[9]};
};

// Per-point functions

inline Tet::Tet(const double* x0, const double* x1, const double* x2, const double* x3)
{
  x0_ = x0[0];
  y0_ = x0[1];
  z0_ = x0[2];

  double j11 = x1[0]-x0[0];
  double j12 = x1[1]-x0[1];
  double j13 = x1[2]-x0[2];
  double j21 = x2[0]-x0[0];
  double j22 = x2[1]-x0[1];
  double j23 = x2[2]-x0[2];
  double j31 = x3[0]-x0[0];
  double j32 = x3[1]-x0[1];
  double j33 = x3[2]-x0[2];

  g2x_ = j22*j33-j23*j32;  g3x_ = j13*j32-j12*j33;  g4x_ = j12*j23-j13*j22;
  g2y_ = j23*j31-j21*j33;  g3y_ = j11*j33-j13*j31;  g4y_ = j13*j21-j11*j23;
  g2z_ = j21*j32-j22*j31;  g3z_ = j12*j31-j11*j32;  g4z_ = j11*j22-j12*j21;
  double det = j11 * g2x_ + j12 * g2y_ + j13 * g2z_;
  double d = 1.0/det;
  g2x_ *= d;  g3x_ *= d;  g4x_ *= d;
  g2y_ *= d;  g3y_ *= d;  g4y_ *= d;
  g2z_ *= d;  g3z_ *= d;  g4z_ *= d;
  g1x_ = -g2x_-g3x_-g4x_;  g1y_ = -g2y_-g3y_-g4y_;  g1z_ = -g2z_-g3z_-g4z_;
}

inline void Tet::xyz2bary(double x, double y, double z,
                          double &r,double &s,double &t,double &u) const
{
  double dx=x-x0_, dy=y-y0_, dz=z-z0_;
  s = g2x_*dx+g2y_*dy+g2z_*dz;
  t = g3x_*dx+g3y_*dy+g3z_*dz;
  u = g4x_*dx+g4y_*dy+g4z_*dz;
  r = 1.-s-t-u;
}

inline bool Tet::contains(const Vector3d& x, std::array<double, 4>& bary) const
{
  xyz2bary(x[0], x[1], x[2], bary[0], bary[1], bary[2], bary[3]);
  return (bary[0] >= 0. && bary[1] >= 0. && bary[2] >= 0. && bary[3] >= 0.);
}

inline Vector3d Tet::bary_grad(const int k) const
{
  switch (k){
    case 0: return {g1x_, g1y_, g1z_};
    case 1: return {g2x_, g2y_, g2z_};
    case 2: return {g3x_, g3y_, g3z_};
    default: return {g4x_, g4y_, g4z_};
  }
}

inline void Tet::linearbasis(double r,
                             double s,
                             double t,
                             double u,
                             double *N) const
{
  N[0] = r;
  N[1] = s;
  N[2] = t;
  N[3] = u;
}

inline void Tet::linearderiv(double ,
                             double ,
                             double ,
                             double ,
                             double *Nx,
                             double *Ny,
                             double *Nz) const
{
  Nx[0] = g1x_;
  Nx[1] = g2x_;
  Nx[2] = g3x_;
  Nx[3] = g4x_;

  Ny[0] = g1y_;
  Ny[1] = g2y_;
  Ny[2] = g3y_;
  Ny[3] = g4y_;

  Nz[0] = g1z_;
  Nz[1] = g2z_;
  Nz[2] = g3z_;
  Nz[3] = g4z_;
}

inline void Tet::quadbasis(double r,
                           double s,
                           double t,
                           double u,
                           double *N) const
{
  N[0] = r*(2*r-1);
  N[1] = s*(2*s-1);
  N[2] = t*(2*t-1);
  N[3] = u*(2*u-1);
  N[perm_[4]] = 4*r*s;
  N[perm_[5]] = 4*s*t;
  N[perm_[6]] = 4*r*t;
  N[perm_[7]] = 4*r*u;
  N[perm_[8]] = 4*s*u;
  N[perm_[9]] = 4*t*u;
}

inline void Tet::quadderiv(double r,
                           double s,
                           double t,
                           double u,
                           double *Nx,
                           double *Ny,
                           double *Nz) const
{
  double a = 4.0*r-1.0;
  double b = 4.0*s-1.0;
  double c = 4.0*t-1.0;
  double d = 4.0*u-1.0;

  Nx[0] = a*g1x_;
  Nx[1] = b*g2x_;
  Nx[2] = c*g3x_;
  Nx[3] = d*g4x_;
  Nx[perm_[4]] = 4*(r*g2x_+s*g1x_);
  Nx[perm_[5]] = 4*(s*g3x_+t*g2x_);
  Nx[perm_[6]] = 4*(t*g1x_+r*g3x_);
  Nx[perm_[7]] = 4*(r*g4x_+u*g1x_);
  Nx[perm_[8]] = 4*(s*g4x_+u*g2x_);
  Nx[perm_[9]] = 4*(t*g4x_+u*g3x_);

  Ny[0] = a*g1y_;
  Ny[1] = b*g2y_;
  Ny[2] = c*g3y_;
  Ny[3] = d*g4y_;
  Ny[perm_[4]] = 4*(r*g2y_+s*g1y_);
  Ny[perm_[5]] = 4*(s*g3y_+t*g2y_);
  Ny[perm_[6]] = 4*(t*g1y_+r*g3y_);
  Ny[perm_[7]] = 4*(r*g4y_+u*g1y_);
  Ny[perm_[8]] = 4*(s*g4y_+u*g2y_);
  Ny[perm_[9]] = 4*(t*g4y_+u*g3y_);

  Nz[0] = a*g1z_;
  Nz[1] = b*g2z_;
  Nz[2] = c*g3z_;
  Nz[3] = d*g4z_;
  Nz[perm_[4]] = 4*(r*g2z_+s*g1z_);
  Nz[perm_[5]] = 4*(s*g3z_+t*g2z_);
  Nz[perm_[6]] = 4*(t*g1z_+r*g3z_);
  Nz[perm_[7]] = 4*(r*g4z_+u*g1z_);
  Nz[perm_[8]] = 4*(s*g4z_+u*g2z_);
  Nz[perm_[9]] = 4*(t*g4z_+u*g3z_);
}

#endif
