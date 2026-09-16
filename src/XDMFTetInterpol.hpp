#ifdef USE_DOLFIN
#ifndef __XDMFTETINTERPOL_HPP
#define __XDMFTETINTERPOL_HPP

#include "Interpol.hpp"
#include "Timestamps.hpp"
#include "Tet.hpp"
#include "cell_locate.hpp"
#include <omp.h>

class XDMFTetInterpol final
  : public Interpol
{
public:
  XDMFTetInterpol(const std::string& infilename);
  ~XDMFTetInterpol() { std::cout << "Destructing XDMFTetInterpol." << std::endl; };
  void update(const double t);
  bool locate(const Vector3d &x, const double t, CellPos& pos);
  void evaluate(const Vector3d &x, const double t, const CellPos& pos, PointValues& );
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
  void print_found() { print_found_counts(found_); }
  Vector3d get_boundary_normal(const Vector3d &x, int & cell_id)
  {
    return cell_normal_[cell_id];
  }
protected:
  MultiTimestamps ts;
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
  bool include_phi = true;

  std::map<std::string, std::string> dolfin_params;

  std::shared_ptr<dolfin::Mesh> mesh;
  Uint dim;

  std::shared_ptr<dolfin::FunctionSpace> u_space_;
  std::shared_ptr<dolfin::FunctionSpace> p_space_;
  // Assuming phi_space_ == p_space_

  // These are new:
  std::vector<double> u_prev_data_;
  std::vector<double> u_next_data_;
  std::vector<double> p_prev_data_;
  std::vector<double> p_next_data_;
  std::vector<double> phi_prev_data_;
  std::vector<double> phi_next_data_;

  std::vector<Tet> tets_;
  std::vector<dolfin::Cell> dolfin_cells_;

  std::vector<CellNeighbours> cell2cells_;
  std::vector<std::uint32_t> dolfin2local_;   // empty: dolfin's cell order


  //std::vector<double> u_prev_coefficients_;
  //std::vector<double> u_next_coefficients_;
  //std::vector<double> p_prev_coefficients_;
  //std::vector<double> p_next_coefficients_;

  //std::vector<double> Nu_, Nux_, Nuy_;
  //std::vector<double> Np_;

  Uint ncoeffs_u;
  Uint ncoeffs_p = 0;   // stays 0 when pressure is ignored

  std::vector<FoundCounts> found_;

  Vector3d _modx(const Vector3d&);

  std::string h5filename_u;
  std::string h5filename_p;
  std::string h5filename_phi;

  std::vector<Uint> i2j;
  std::vector<Uint> j2i;

  CellDofs u_dofs_;
  CellDofs p_dofs_;

  std::vector<int> cell_type_;
  std::vector<Vector3d> cell_normal_;
  std::vector<Vector3d> cell_facet_midpoint_;
  std::vector<std::vector<int>> perm_;

};

#endif
#endif
