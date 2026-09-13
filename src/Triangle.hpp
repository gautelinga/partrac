#ifdef USE_DOLFIN
#ifndef __TRIANGLE_HPP
#define __TRIANGLE_HPP

#include <dolfin.h>
#include <array>
#include "typedefs.hpp"

class Triangle
{

public:

  static constexpr std::size_t n_dofs_max = 6;   // quadbasis writes this many

  Triangle() {}
  Triangle(const dolfin::Cell& cell);

  void xy2bary(double x, double y,
               double &r, double &s, double &t) const;

  bool contains(const Vector3d& x) const;

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

  double get_det() const { return det; };
  void dump();
  std::vector<double> dof_coords(const int index) const;
  double dot_grad_gi(const double vx, const double vy, const int index) const;

private:

  std::array<double, 3> xx_, yy_;
  double g1x_, g1y_;
  double g2x_, g2y_;
  double g3x_, g3y_;
  double det;

  static constexpr std::array<int, 6> perm_ = {-1, -1, -1, 5, 3, 4};  // Check!
  // static constexpr std::array<int, 6> perm_alt_ = {-1, -1, -1, 4, 5, 3};
};

#endif
#endif
