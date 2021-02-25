#ifdef USE_DOLFIN
#ifndef TRIANGLE_HPP
#define TRIANGLE_HPP

#include <dolfin.h>
#include <array>

class Triangle
{

public:

  Triangle() {}
  Triangle(const dolfin::Cell& cell);

  void xy2bary(double x, double y,
               double &r, double &s, double &t) const;

  void linearbasis(double r, double s, double t,
                   std::array<double, 3> &N) const;

  void linearderiv(double r, double s, double t,
                   std::array<double, 3> &Nx,
                   std::array<double, 3> &Ny) const;

  void quadbasis(double r, double s, double t,
                 std::array<double, 6> &N) const;

  void quadderiv(double r, double s, double t,
                 std::array<double, 6> &Nx,
                 std::array<double, 6> &Ny) const;

private:

  std::array<double, 3> xx_, yy_;
  double g1x_,g1y_;
  double g2x_,g2y_;
  double g3x_,g3y_;

  static constexpr std::array<int, 6> perm_ = {-1, -1, -1, 5, 3, 4};  // Check!
};

#endif
#endif
