#ifdef USE_DOLFIN
#ifndef __TRIANGLEINTERPOL_HPP
#define __TRIANGLEINTERPOL_HPP

#include "Interpol.hpp"
#include "Timestamps.hpp"
#include "Triangle.hpp"
#include <numeric>
#include <omp.h>

class TriangleInterpol
  : public Interpol
{
public:
  TriangleInterpol(const std::string& infilename);
  ~TriangleInterpol() { std::cout << "Destructing TriangleInterpol." << std::endl; };
  void update(const double t);
  bool locate(const Vector3d &x, const double t, int& id_prev);
  void evaluate(const Vector3d &x, const double t, const int id, PointValues& );
  double get_t_min() { return ts.get_t_min(); };
  double get_t_max() { return ts.get_t_max(); };
  double get_rho() {
    if (contains(dolfin_params, std::string("rho")))
      return stod(dolfin_params["rho"]);
    else {
      std::cout << "dolfin_params does not contain \"rho\"" << std::endl;
      exit(1);
    }
  };
  using Interpol::locate;
  using Interpol::evaluate;
  void print_found() {
    auto found_same = std::reduce(found_same_.begin(), found_same_.end());
    auto found_nneigh = std::reduce(found_nneigh_.begin(), found_nneigh_.end());
    auto found_other = std::reduce(found_other_.begin(), found_other_.end());

    auto found_sum = found_same + found_nneigh + found_other;
    double frac_same = double(found_same) / found_sum;
    double frac_nneigh = double(found_nneigh) / found_sum;
    double frac_other = 1. - frac_same - frac_nneigh;
    std::cout << "Found in same cell: " << frac_same << ", nearest neighbour cell: " << frac_nneigh << ", other cell: " << frac_other << std::endl;
    std::fill(found_same_.begin(), found_same_.end(), 0);
    std::fill(found_nneigh_.begin(), found_nneigh_.end(), 0);
    std::fill(found_other_.begin(), found_other_.end(), 0);
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

  std::shared_ptr<dolfin::Function> u_prev_;
  std::shared_ptr<dolfin::Function> u_next_;
  std::shared_ptr<dolfin::Function> p_prev_;
  std::shared_ptr<dolfin::Function> p_next_;

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

  // One counter per thread, to avoid races in locate
  std::vector<long unsigned int> found_same_;
  std::vector<long unsigned int> found_nneigh_;
  std::vector<long unsigned int> found_other_;

  void _modx(dolfin::Array<double>&, const Vector3d&);

};

#endif
#endif
