#ifndef __POINTVALUES_HPP
#define __POINTVALUES_HPP

#include <array>
#include "typedefs.hpp"

// A point's cell and its barycentric coordinates there
struct CellPos {
  int id = -1;
  std::array<double, 4> bary;
};

class PointValues {
public:
  PointValues(const double U0) : U0(U0) {
    U = {0., 0., 0.};
    A = {0., 0., 0.};
    gradU << 0., 0., 0., 0., 0., 0., 0., 0., 0.; 
    gradA << 0., 0., 0., 0., 0., 0., 0., 0., 0.; 
    P = 0.;
    Rho = 0.;
  };
  Vector3d U;
  Vector3d A;
  Matrix3d gradU;
  Matrix3d gradA;
  double P = 0.;
  double Rho = 0.;    // density
  double Phi = 0.;    // phase field (XDMF)
  Vector3d get_u() { return U0 * U; };
  Matrix3d get_J() { return U0 * gradU; };
  Vector3d get_Ju() { return U0 * U0 * gradU * U; }; // check
  Vector3d get_a() { return U0 * A; };
  Matrix3d get_grada() { return U0 * gradA; }
  double get_p() const { return P; };
  double get_rho() const { return Rho; };
  double get_phi() const { return Phi; };
  int cell_type = 0;
  int get_cell_type() const { return cell_type; };
private:
  double U0;
};

#endif
