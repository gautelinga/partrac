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
#include "geometry.hpp"
#include "strings.hpp"

template<typename Cell>
class MeshInterpol : public Interpol {
public:
  MeshInterpol(const std::string& infilename) : Interpol(infilename) {}
  bool locate(const Vector3d &x, const double t, CellPos& pos){
    assert(t <= t_next && t >= t_prev);
    const Vector3d xx = _modx(x);
    return walk_to_cell(cells_, facet_neigh_, xx, pos) || locate_tree(xx, pos);
  }
  bool reflect(const Vector3d& x, Vector3d& dx, CellPos& pos);
  void enable_reflection();
  // Outward unit normal of the cell's wall facets, their mean at an edge; zero off the wall
  Vector3d get_boundary_normal(const Vector3d &x, int& cell_id);
  double hmin() const { return mesh->hmin(); }
  using Interpol::locate;
  using Interpol::evaluate;
protected:
  // Periodicity and the pressure flag from the parameter file
  void read_mesh_params();
  // Dimension, mesh tables, bounding box tree and the domain bounds
  void init_mesh_geometry();
  // Cells in dof order, with their facet table
  void build_cells(const dolfin::GenericDofMap& dofmap);
  // Out of the step loops: the tree, when the walk does not find the cell
  bool locate_tree(const Vector3d& xx, CellPos& pos);
  // Into the box along the periodic axes
  Vector3d _modx(const Vector3d& x) const {
    Vector3d x_loc = x;
    for (Uint i = 0; i < dim; ++i)
      if (period_[i] > 0. && (x[i] < x_min[i] || x[i] >= x_max[i]))
        x_loc[i] = x_min[i] + modulox(x[i] - x_min[i], period_[i]);
    return x_loc;
  }
  // The neighbour across each facet: a cell, a wall or a periodic image
  void build_facet_table();

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
  std::vector<std::uint32_t> dolfin2local_;   // empty: dolfin's cell order
  std::vector<std::int32_t> facet_neigh_;     // walk_to_cell, reflect_in_cells
  Vector3d period_ = Vector3d::Zero();
  CellDofs u_dofs_;
  CellDofs p_dofs_;
};

#endif
#endif
