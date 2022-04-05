#ifdef USE_DOLFIN
#ifndef __TRIANGLEINTERPOL_HPP
#define __TRIANGLEINTERPOL_HPP

#include "Interpol.hpp"
#include "Timestamps.hpp"
#include "Triangle.hpp"

class TriangleInterpol
  : public Interpol
{
public:
  TriangleInterpol(const std::string& infilename);
  void update(const double t);
  void probe(const Vector3d &x, const double t);
  bool inside_domain() const { return inside; };
  double get_ux(){ return U[0]; };
  double get_uy(){ return U[1]; };
  double get_uz(){ return 0.; };
  double get_ax() { return A[0]; };
  double get_ay() { return A[1];};
  double get_az() { return 0.;};
  double get_t_min() { return ts.get_t_min(); };
  double get_t_max() { return ts.get_t_max(); };
  double get_rho() {
    if (contains(dolfin_params, std::string("rho")))
      return stod(dolfin_params["rho"]);
    else {
      std::cout << "dolfin_params does not contain \"rho\"" << std::endl;
      exit(0);
    }
  };
  double get_p() { return P;};
  double get_uxx() { return gradU(0, 0); };
  double get_uxy() { return gradU(0, 1); };
  double get_uxz() { return 0.; };
  double get_uyx() { return gradU(1, 0); };
  double get_uyy() { return gradU(1, 1); };
  double get_uyz() { return 0.; };
  double get_uzx() { return 0.; };
  double get_uzy() { return 0.; };
  double get_uzz() { return 0.; };
  Matrix3d get_grada() { return gradA; };

protected:
  Timestamps ts;
  double t_prev = 0.;
  double t_next = 0.;
  double alpha_t;

  std::vector<bool> periodic = {false, false, false};
  //Vector3d x_min = {0., 0., 0.};
  //Vector3d x_max = {0., 0., 0.};

  Vector3d U = {0., 0., 0.};  // FIXME: 2d
  //double Uy = 0.;
  //double Uz = 0.;
  Vector3d A = {0., 0., 0.};
  double P = 0.;
  Matrix3d gradU, gradA;

  bool inside;

  bool include_pressure = true;

  std::map<std::string, std::string> dolfin_params;

  std::shared_ptr<dolfin::Mesh> mesh;
  Uint dim;

  std::shared_ptr<dolfin::FunctionSpace> u_space_;
  std::shared_ptr<dolfin::FunctionSpace> p_space_;

  std::shared_ptr<dolfin::Function> u_prev_;
  std::shared_ptr<dolfin::Function> u_next_;
  std::shared_ptr<dolfin::Function> p_prev_;
  std::shared_ptr<dolfin::Function> p_next_;

  std::vector<Triangle> triangles_;
  std::vector<dolfin::Cell> dolfin_cells_;
  std::vector<ufc::cell> ufc_cells_;

  std::vector<std::vector<double>> coordinate_dofs_;

  //std::array<double, 12> u_prev_coefficients_;
  //std::array<double, 12> u_next_coefficients_;
  //std::array<double, 3> p_prev_coefficients_;
  //std::array<double, 3> p_next_coefficients_;
  std::vector<double> u_prev_coefficients_;
  std::vector<double> u_next_coefficients_;
  std::vector<double> p_prev_coefficients_;
  std::vector<double> p_next_coefficients_;

  //std::array<double, 6> N6_, Nx_, Ny_, Nz_;
  //std::array<double, 3> N3_;
  std::vector<double> Nu_, Nux_, Nuy_;
  std::vector<double> Np_;

  Uint ncoeffs_u;
  Uint ncoeffs_p;

};

#endif
#endif
