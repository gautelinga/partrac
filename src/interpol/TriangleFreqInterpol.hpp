#ifdef USE_DOLFIN
#ifndef __TRIANGLEFREQINTERPOL_HPP
#define __TRIANGLEFREQINTERPOL_HPP

#include "Interpol.hpp"
#include "Params.hpp"
#include "strings.hpp"
#include "FreqStamps.hpp"
#include "Triangle.hpp"
#include "cell_locate.hpp"

class TriangleFreqInterpol final
  : public Interpol
{
public:
  TriangleFreqInterpol(const std::string& infilename);
  ~TriangleFreqInterpol() { std::cout << "Destructing TriangleFreqInterpol." << std::endl; };
  void update(const double t);
  bool locate(const Vector3d &x, const double t, CellPos& pos);
  void evaluate(const Vector3d &x, const double t, const CellPos& pos, PointValues& fields);
  double get_t_min() { return dolfin_params.get<double>("t_min"); };
  double get_t_max() { return dolfin_params.get<double>("t_max"); };
  double get_rho() {
    if (dolfin_params.has("rho"))
      return dolfin_params.get<double>("rho");
    else {
      std::cout << "dolfin_params does not contain \"rho\"" << std::endl;
      exit(1);
    }
  };
  using Interpol::locate;
  using Interpol::evaluate;
  // Walk a move off the walls
  bool reflect(const Vector3d& x, Vector3d& dx, CellPos& pos);
  void enable_reflection();
  double hmin() const { return mesh->hmin(); }
  void print_found() { print_found_counts(found_); }
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

  partrac::Params dolfin_params;

  std::shared_ptr<dolfin::Mesh> mesh;
  Uint dim;

  std::shared_ptr<dolfin::FunctionSpace> u_space_;
  std::shared_ptr<dolfin::FunctionSpace> p_space_;

  std::vector<Triangle> triangles_;
  std::vector<dolfin::Cell> dolfin_cells_;

  std::vector<CellNeighbours> cell2cells_;
  std::vector<std::int32_t> facet_neigh_;   // reflect_in_cells
  Vector3d period_ = Vector3d::Zero();
  std::vector<std::uint32_t> dolfin2local_;   // empty: dolfin's cell order

  //std::vector<double> Nu_, Nux_, Nuy_;
  //std::vector<double> Np_;

  Uint ncoeffs_u;
  Uint ncoeffs_p = 0;   // stays 0 when pressure is ignored

  std::shared_ptr<dolfin::Function> u_;
  std::shared_ptr<dolfin::Function> p_;

  /*
  std::vector<std::shared_ptr<dolfin::Function>> u__;
  std::vector<std::shared_ptr<dolfin::Function>> p__;
  std::vector<std::vector<double>> u_coefficients__;
  std::vector<std::vector<double>> p_coefficients__;
  */ 
 
  std::vector<std::vector<std::vector<double>>> u_coefficients_;
  std::vector<std::vector<std::vector<double>>> p_coefficients_;

  std::vector<FoundCounts> found_;

  double omega0 = 0.;

  Vector3d _modx(const Vector3d&);

};

#endif
#endif
