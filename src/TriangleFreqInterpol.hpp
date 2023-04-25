#ifdef USE_DOLFIN
#ifndef __TRIANGLEFREQINTERPOL_HPP
#define __TRIANGLEFREQINTERPOL_HPP

#include "Interpol.hpp"
#include "FreqStamps.hpp"
#include "Triangle.hpp"

class TriangleFreqInterpol
  : public Interpol
{
public:
  TriangleFreqInterpol(const std::string& infilename);
  ~TriangleFreqInterpol() { std::cout << "Destructing TriangleFreqInterpol." << std::endl; };
  void update(const double t);
  void probe(const Vector3d &x, const double t);
  void probe(const Vector3d &x, const double t, int& id_prev);
  bool inside_domain() const { return inside; };
  double get_ux(){ return U[0]; };
  double get_uy(){ return U[1]; };
  double get_uz(){ return 0.; };
  double get_ax() { return A[0]; };
  double get_ay() { return A[1];};
  double get_az() { return 0.;};
  double get_t_min() { return stod(dolfin_params["t_min"]); };
  double get_t_max() { return stod(dolfin_params["t_max"]); };
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
  void probe(const Vector3d &x){ probe(x, this->t_update); };
  void print_found() {
    long int found_sum = found_same + found_nneigh + found_other;
    double frac_same = double(found_same) / found_sum;
    double frac_nneigh = double(found_nneigh) / found_sum;
    double frac_other = 1. - frac_same - frac_nneigh;
    std::cout << "Found in same cell: " << frac_same << ", nearest neighbour cell: " << frac_nneigh << ", other cell: " << frac_other << std::endl;
    found_same = 0;
    found_nneigh = 0;
    found_other = 0;
  }
protected:
  FreqStamps fs; // frequencies holder
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

  std::vector<Triangle> triangles_;
  std::vector<dolfin::Cell> dolfin_cells_;
  std::vector<ufc::cell> ufc_cells_;

  std::vector<std::set<Uint>> cell2cells_;
  std::vector<std::vector<double>> coordinate_dofs_;

  std::vector<double> Nu_, Nux_, Nuy_;
  std::vector<double> Np_;

  Uint ncoeffs_u;
  Uint ncoeffs_p;

  std::vector<std::shared_ptr<dolfin::Function>> u__;
  std::vector<std::shared_ptr<dolfin::Function>> p__;

  //std::shared_ptr<dolfin::Function> p_1_;
  //std::shared_ptr<dolfin::Function> u_1_;

  std::vector<std::vector<double>> u_coefficients__;
  std::vector<std::vector<double>> p_coefficients__;

  //std::vector<double> p_1_coefficients_;
  //std::vector<double> u_1_coefficients_;

  long unsigned int found_same = 0;
  long unsigned int found_nneigh = 0;
  long unsigned int found_other = 0;

  std::vector<double> w_f_;
  std::vector<double> wt_f_;

  std::vector<double> ux_f_;
  std::vector<double> uy_f_;
  std::vector<double> p_f_;

  std::vector<double> uxx_f_;
  std::vector<double> uxy_f_;
  std::vector<double> uyx_f_;
  std::vector<double> uyy_f_;

  double omega0 = 0.;
};

#endif
#endif
