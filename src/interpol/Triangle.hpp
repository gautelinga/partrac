#ifndef __TRIANGLE_HPP
#define __TRIANGLE_HPP

#include <array>
#include "typedefs.hpp"

#ifdef USE_DOLFIN
namespace dolfin { class Cell; }
#endif

class alignas(64) Triangle
{

public:

  static constexpr int n_verts = 3;
  static constexpr std::size_t n_dofs_max = 6;   // quadbasis writes this many

  Triangle() {}
  // The three vertex coordinate pairs, in dolfin's local vertex order
  Triangle(const double* x0, const double* x1, const double* x2);
#ifdef USE_DOLFIN
  Triangle(const dolfin::Cell& cell);
#endif

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

inline Triangle::Triangle(const double* x0, const double* x1, const double* x2)
{
  x0_ = x0[0];
  y0_ = x0[1];

  double j11 = x1[0]-x0[0];
  double j12 = x1[1]-x0[1];
  double j21 = x2[0]-x0[0];
  double j22 = x2[1]-x0[1];

  const double det = j11 * j22 - j12*j21;

  double d = 1.0/det;
  g2x_ = j22*d;   g3x_ = -j12*d;
  g2y_ = -j21*d;  g3y_ = j11*d;
  g1x_ = -g2x_-g3x_;  g1y_ = -g2y_-g3y_;
}

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
