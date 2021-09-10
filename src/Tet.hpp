#ifdef USE_DOLFIN
#ifndef TET_HPP
#define TET_HPP

#include <dolfin.h>
#include <array>

class Tet
{

public:

  Tet() {}
  Tet(const dolfin::Cell& cell);

  void xyz2bary(double x, double y, double z,
                double &r, double &s, double &t, double &u) const;

  void linearbasis(double r, double s, double t, double u,
                   //std::array<double, 4> &N) const;
                   std::vector<double> &N) const;

  void linearderiv(double, double, double, double,
                   //std::array<double, 4> &Nx,
                   //std::array<double, 4> &Ny,
                   //std::array<double, 4> &Nz) const;
                   std::vector<double> &Nx,
                   std::vector<double> &Ny,
                   std::vector<double> &Nz) const;

  void quadbasis(double r, double s, double t, double u,
                 //std::array<double, 10> &N) const;
                 std::vector<double> &N) const;

  void quadderiv(double r,double s,double t,double u,
                 //std::array<double, 10> &Nx,
                 //std::array<double, 10> &Ny,
                 //std::array<double, 10> &Nz) const;
                 std::vector<double> &Nx,
                 std::vector<double> &Ny,
                 std::vector<double> &Nz) const;

private:

  std::array<double, 4> xx_, yy_, zz_;
  double g1x_, g1y_, g1z_;
  double g2x_, g2y_, g2z_;
  double g3x_, g3y_, g3z_;
  double g4x_, g4y_, g4z_;

  static constexpr std::array<int, 10> perm_ = {-1, -1, -1, -1, 9, 6, 8, 7, 5, 4};
};

#endif
#endif
