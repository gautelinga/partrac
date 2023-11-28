#ifdef USE_DOLFIN
#ifndef __XDMFTRIANGLEINTERPOL_HPP
#define __XDMFTRIANGLEINTERPOL_HPP

#include "Interpol.hpp"
#include "Timestamps.hpp"
#include "Triangle.hpp"
#include <omp.h>

class XDMFTriangleInterpol
  : public Interpol
{
public:
  XDMFTriangleInterpol(const std::string& infilename);
  ~XDMFTriangleInterpol() { std::cout << "Destructing XDMFTriangleInterpol." << std::endl; };
  void update(const double t);
  void probe(const Vector3d &x, const double t);
  void probe(const Vector3d &x, const double t, int& id_prev);
  bool probe_light(const Vector3d &x, const double t, int& id_prev);
  void probe_heavy(const Vector3d &x, const double t, const int id, PointValues& );
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
  void probe(const Vector3d &x){ probe(x, this->t_update); };
  void print_found()
  {
    auto found_same = std::reduce(found_same_.begin(), found_same_.end());
    auto found_nneigh = std::reduce(found_nneigh_.begin(), found_nneigh_.end());
    auto found_other = std::reduce(found_other_.begin(), found_other_.end());

    auto found_sum = found_same + found_nneigh + found_other;
    auto frac_same = double(found_same) / found_sum;
    auto frac_nneigh = double(found_nneigh) / found_sum;
    auto frac_other = 1. - frac_same - frac_nneigh;
    std::cout << "Found in same cell: " << frac_same << ", nearest neighbour cell: " << frac_nneigh << ", other cell: " << frac_other << std::endl;
    found_same = 0;
    found_nneigh = 0;
    found_other = 0;
  }
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

  // These are new:
  std::vector<double> u_prev_data_;
  std::vector<double> u_next_data_;
  std::vector<double> p_prev_data_;
  std::vector<double> p_next_data_;

  // End

  std::vector<Triangle> triangles_;
  std::vector<dolfin::Cell> dolfin_cells_;
  std::vector<ufc::cell> ufc_cells_;

  std::vector<std::set<Uint>> cell2cells_;

  std::vector<std::vector<double>> coordinate_dofs_;

  std::vector<double> u_prev_coefficients_;
  std::vector<double> u_next_coefficients_;
  std::vector<double> p_prev_coefficients_;
  std::vector<double> p_next_coefficients_;

  std::vector<double> Nu_, Nux_, Nuy_;
  std::vector<double> Np_;

  Uint ncoeffs_u;
  Uint ncoeffs_p;

  std::vector<long unsigned int> found_same_;
  std::vector<long unsigned int> found_nneigh_;
  std::vector<long unsigned int> found_other_;

  void _modx(dolfin::Array<double>&, const Vector3d&);

  std::string h5filename_u;
  std::string h5filename_p;
  std::string h5filename_phi;

  std::vector<Uint> i2j;
  std::vector<Uint> j2i;

};

#endif
#endif
