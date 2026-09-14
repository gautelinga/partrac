#ifdef USE_DOLFIN
#ifndef __TETINTERPOL_HPP
#define __TETINTERPOL_HPP

#include "Tet.hpp"
#include "cell_locate.hpp"
#include "Interpol.hpp"
#include "Timestamps.hpp"

class TetInterpol final
  : public Interpol
{
public:

  TetInterpol(const std::string& infilename);
  void update(const double t);
  using Interpol::locate;
  using Interpol::evaluate;
  bool locate(const Vector3d &x, const double t, int& cell_id);
  void evaluate(const Vector3d &x, const double t, const int cell_id, PointValues& ptvals);
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
  void print_found() { print_found_counts(found_); }
  void reflect(Vector3d &x, Vector3d &dx_new, const double t, const double dt, int& cell_id);
protected:

  Timestamps ts;
  double t_prev = 0.;
  double t_next = 0.;
  // double alpha_t;

  std::vector<bool> periodic = {false, false, false};
  //Vector3d x_min = {0., 0., 0.};
  //Vector3d x_max = {0., 0., 0.};

  Vector3d U = {0., 0., 0.};
  //double Uy = 0.;
  //double Uz = 0.;
  Vector3d A = {0., 0., 0.};
  double P = 0.;
  Matrix3d gradU, gradA;

  bool inside;

  std::map<std::string, std::string> dolfin_params;

  bool include_pressure = true;

  std::shared_ptr<dolfin::Mesh> mesh;
  Uint dim;

  std::shared_ptr<dolfin::FunctionSpace> u_space_;
  std::shared_ptr<dolfin::FunctionSpace> p_space_;

  std::shared_ptr<dolfin::Function> u_prev_;
  std::shared_ptr<dolfin::Function> u_next_;
  std::shared_ptr<dolfin::Function> p_prev_;
  std::shared_ptr<dolfin::Function> p_next_;

  std::vector<double> u_prev_vec;
  std::vector<double> u_next_vec;
  std::vector<double> p_prev_vec;
  std::vector<double> p_next_vec;

  std::vector<Tet> tets_;
  std::vector<dolfin::Cell> dolfin_cells_;


  std::vector<CellNeighbours> cell2cells_;

  // std::array<double, 30> u_prev_coefficients_;
  // std::array<double, 30> u_next_coefficients_;
  // std::array<double, 4> p_prev_coefficients_;
  // std::array<double, 4> p_next_coefficients_;
  
  // std::vector<double> u_prev_coefficients_;
  // std::vector<double> u_next_coefficients_;
  // std::vector<double> p_prev_coefficients_;
  // std::vector<double> p_next_coefficients_;

  // std::array<double, 10> N10_, Nx_, Ny_, Nz_;
  // std::array<double, 4> N4_;
  // std::vector<double> Nu_, Nux_, Nuy_, Nuz_;
  // std::vector<double> Np_;

  Uint ncoeffs_u;
  Uint ncoeffs_p = 0;   // stays 0 when pressure is ignored

  std::vector<FoundCounts> found_;

  CellDofs u_dofs_;
  CellDofs p_dofs_;

  Vector3d _modx(const Vector3d&);

  std::vector<int> cell_type_;
  std::vector<std::vector<int>> cell_facets_;
  std::vector<std::vector<Vector3d>> facets_;

  double hmin;
  bool _cross_facet(double& beta, Vector3d& N, const Vector3d& x, const Vector3d &dx_new, std::vector<Vector3d> &facet);
};

#endif
#endif
