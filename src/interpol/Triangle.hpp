#ifdef USE_DOLFIN
#ifndef __TRIANGLE_HPP
#define __TRIANGLE_HPP

#include <dolfin.h>
#include <array>
#include "typedefs.hpp"

class alignas(64) Triangle
{

public:

  static constexpr int n_verts = 3;
  static constexpr std::size_t n_dofs_max = 6;   // quadbasis writes this many

  Triangle() {}
  Triangle(const dolfin::Cell& cell);

  void xy2bary(double x, double y,
               double &r, double &s, double &t) const;

  // Writes the first three barycentrics, inside or not
  bool contains(const Vector3d& x, std::array<double, 4>& bary) const;

  // Gradient of barycentric k: inward normal of the edge opposite vertex k
  Vector3d bary_grad(const int k) const;

  void linearbasis( double r
                  , double s
                  , double t
                  , double *N
                  ) const;

  void linearderiv( double r
                  , double s
                  , double t
                  , double *Nx
                  , double *Ny
                  ) const;

  void quadbasis( double r
                , double s
                , double t
                , double *N
                ) const;

  void quadderiv( double r
                , double s
                , double t
                , double *Nx
                , double *Ny
                ) const;


private:

  // Basis data only
  double x0_, y0_;
  double g2x_, g2y_;
  double g3x_, g3y_;
  double g1x_, g1y_;

  static constexpr std::array<int, 6> perm_ = {-1, -1, -1, 5, 3, 4};  // Check!
  // static constexpr std::array<int, 6> perm_alt_ = {-1, -1, -1, 4, 5, 3};

public:

  // quadbasis slots of the midpoints of edges 01, 02, 12
  static constexpr std::array<int, 3> mid_ = {perm_[3], perm_[5], perm_[4]};
};

// Per-point functions

inline void Triangle::xy2bary(double x, double y,
                              double &r, double &s, double &t) const
{
  double dx=x-x0_, dy=y-y0_;
  s = g2x_*dx+g2y_*dy;
  t = g3x_*dx+g3y_*dy;
  r = 1.-s-t;
}

inline bool Triangle::contains(const Vector3d& x, std::array<double, 4>& bary) const
{
  xy2bary(x[0], x[1], bary[0], bary[1], bary[2]);
  return (bary[0] >= 0. && bary[1] >= 0. && bary[2] >= 0.);
}

inline Vector3d Triangle::bary_grad(const int k) const
{
  switch (k){
    case 0: return {g1x_, g1y_, 0.};
    case 1: return {g2x_, g2y_, 0.};
    default: return {g3x_, g3y_, 0.};
  }
}

inline void Triangle::linearbasis( double r
                                 , double s
                                 , double t
                                 , double *N
                                 ) const
{
  N[0] = r;
  N[1] = s;
  N[2] = t;
}

inline void Triangle::linearderiv( double r
                                 , double s
                                 , double t
                                 , double *Nx
                                 , double *Ny
                                 ) const
{
  Nx[0] = g1x_;
  Nx[1] = g2x_;
  Nx[2] = g3x_;

  Ny[0] = g1y_;
  Ny[1] = g2y_;
  Ny[2] = g3y_;
}

inline void Triangle::quadbasis( double r
                               , double s
                               , double t
                               , double *N
                               ) const
{
  N[0] = r*(2*r-1);
  N[1] = s*(2*s-1);
  N[2] = t*(2*t-1);
  N[perm_[3]] = 4*r*s;
  N[perm_[4]] = 4*s*t;
  N[perm_[5]] = 4*r*t;
}

inline void Triangle::quadderiv( double r
                               , double s
                               , double t
                               , double *Nx
                               , double *Ny
                               ) const
{
  double a = 4.0*r-1.0;
  double b = 4.0*s-1.0;
  double c = 4.0*t-1.0;

  Nx[0] = a*g1x_;
  Nx[1] = b*g2x_;
  Nx[2] = c*g3x_;
  Nx[perm_[3]] = 4*(r*g2x_+s*g1x_);
  Nx[perm_[4]] = 4*(s*g3x_+t*g2x_);
  Nx[perm_[5]] = 4*(t*g1x_+r*g3x_);

  Ny[0] = a*g1y_;
  Ny[1] = b*g2y_;
  Ny[2] = c*g3y_;
  Ny[perm_[3]] = 4*(r*g2y_+s*g1y_);
  Ny[perm_[4]] = 4*(s*g3y_+t*g2y_);
  Ny[perm_[5]] = 4*(t*g1y_+r*g3y_);
}

#endif
#endif
