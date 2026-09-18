#ifdef USE_DOLFIN
#ifndef __MESHINTERPOL_HPP
#define __MESHINTERPOL_HPP

// What every interpolator on a simplex mesh holds and does: the cells and
// their tables, locating a point, and walking a move off the walls. The
// loaders add the fields: their own update and evaluate.

#include <memory>
#include <string>
#include <vector>
#include <dolfin.h>
#include "Interpol.hpp"
#include "Params.hpp"
#include "cell_locate.hpp"
#include "strings.hpp"

template<typename Cell>
class MeshInterpol : public Interpol {
public:
  MeshInterpol(const std::string& infilename) : Interpol(infilename) {}
  bool locate(const Vector3d &x, const double t, CellPos& pos);
  bool reflect(const Vector3d& x, Vector3d& dx, CellPos& pos);
  void enable_reflection();
  double hmin() const { return mesh->hmin(); }
  void print_found() { print_found_counts(found_); }
  double get_rho() {
    if (dolfin_params.has("rho"))
      return dolfin_params.get<double>("rho");
    std::cout << "dolfin_params does not contain \"rho\"" << std::endl;
    exit(1);
  };
  using Interpol::locate;
  using Interpol::evaluate;
protected:
  // Periodicity and the pressure flag from the parameter file
  void read_mesh_params();
  // Dimension, mesh tables, bounding box tree and the domain bounds
  void init_mesh_geometry();
  // Cells in dof order, with their neighbours; pair_periodic for periodic neighbours
  void build_cells(const dolfin::GenericDofMap& dofmap, const bool pair_periodic);
  Vector3d _modx(const Vector3d&);

  partrac::Params dolfin_params;
  double t_prev = 0.;   // the stamps the fields are between
  double t_next = 0.;
  std::vector<bool> periodic = {false, false, false};
  bool include_pressure = true;
  double periodic_tol = 1e-12;   // heuristic

  std::shared_ptr<dolfin::Mesh> mesh;
  Uint dim;
  std::shared_ptr<dolfin::FunctionSpace> u_space_;
  std::shared_ptr<dolfin::FunctionSpace> p_space_;
  Uint ncoeffs_u;
  Uint ncoeffs_p = 0;   // stays 0 when pressure is ignored

  std::vector<Cell> cells_;
  std::vector<dolfin::Cell> dolfin_cells_;
  std::vector<CellNeighbours> cell2cells_;
  std::vector<std::uint32_t> dolfin2local_;   // empty: dolfin's cell order
  std::vector<std::int32_t> facet_neigh_;     // reflect_in_cells
  Vector3d period_ = Vector3d::Zero();
  CellDofs u_dofs_;
  CellDofs p_dofs_;
  std::vector<FoundCounts> found_;
};

#endif
#endif
